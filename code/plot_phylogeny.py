#!/usr/bin/env python3
"""Render an IQ-TREE divergence phylogeny with Baltic 1.x and midpoint rooting."""
import argparse
import csv
from importlib.metadata import version
import json
import math
from pathlib import Path
from urllib.parse import unquote
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import baltic as bt
from baltic.bt_utils import plot_scale_bar
from sintax_parser import taxon_parts
from sequence_headers import local_header_parts
from genbank_xml import iter_genbank_records


def split_keys(tree):
    """Canonical unrooted bipartitions, so supports survive changes of root."""
    descendants = {}
    stack = [(tree.root, False)]
    while stack:
        branch, visited = stack.pop()
        if branch.is_leaf():
            descendants[branch] = frozenset([branch.name])
        elif visited:
            descendants[branch] = frozenset().union(*(descendants[c] for c in branch.children))
        else:
            stack.append((branch, True))
            stack.extend((c, False) for c in branch.children)
    all_tips = descendants[tree.root]
    keys = {}
    for branch, side in descendants.items():
        if branch is tree.root:
            continue
        other = all_tips - side
        keys[branch] = min(tuple(sorted(side)), tuple(sorted(other)), key=lambda x: (len(x), x))
    return keys


def midpoint_with_support(tree):
    supports = {}
    for branch, key in split_keys(tree).items():
        label = str(branch.traits.get('label', ''))
        if not branch.is_leaf() and label:
            if key in supports and supports[key] != label:
                raise ValueError(f'Conflicting support labels for the same split: {supports[key]}, {label}')
            supports[key] = label
    lengths = [b.length for b in tree.Objects if b is not tree.root]
    if any(not math.isfinite(x) or x < 0 for x in lengths) or not any(x > 0 for x in lengths):
        raise ValueError('Midpoint rooting requires finite nonnegative branch lengths and a positive tree span')
    tree.midpoint_root()
    for branch in tree.Objects:
        branch.traits.pop('label', None)
    assigned = set()
    for branch, key in split_keys(tree).items():
        if key in supports and key not in assigned:
            branch.traits['label'] = supports[key]
            assigned.add(key)
    if assigned != set(supports):
        raise ValueError('Rerooting lost a supported bipartition')
    return supports


def reference_names(path):
    names = {}
    if Path(path).suffix == '.tsv':
        with Path(path).open() as handle:
            return {row['sequence_id']: row['taxon'] for row in csv.DictReader(handle, delimiter='\t')}
    for record in iter_genbank_records(path):
        organism = record.findtext('GBSeq_organism')
        accession = record.findtext('GBSeq_accession-version')
        primary = record.findtext('GBSeq_primary-accession')
        if organism:
            for identifier in (accession, primary):
                if identifier:
                    names[identifier] = organism
    return names



# Resolve cluster seed labels from analysis tables or external-source metadata.
def tip_metadata(tip_ids, analysis_dir, reference_genbank):
    references = reference_names(reference_genbank)
    bold_codes = {}
    if Path(reference_genbank).suffix == '.tsv':
        with Path(reference_genbank).open() as handle:
            bold_codes = {row['sequence_id']: row['process_id'] for row in csv.DictReader(handle, delimiter='\t')}
    metadata, local = {}, {}
    for identifier in tip_ids:
        parts = local_header_parts(identifier)
        if parts is not None:
            folder, readid = parts
            folder = unquote(folder)
            local.setdefault(folder, {})[readid] = identifier
        else:
            # Modern efetch headers use accession.version; also handle legacy pipe IDs.
            candidates = [identifier] + identifier.split('|')
            accession = next((x for x in candidates if x in references), None)
            if accession is None:
                raise ValueError(f'No source metadata for reference {identifier}')
            organism = references[accession]
            display_id = bold_codes.get(accession, accession)
            metadata[identifier] = [identifier, 'bold' if bold_codes else 'reference', organism, display_id, '', '', f'{organism} | {display_id}']
    for folder, wanted in local.items():
        path = Path(analysis_dir) / folder / f'{folder}_reads.tsv'
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                identifier = wanted.pop(row['readid'], None)
                if identifier is None:
                    continue
                species = taxon_parts(row['accepted_species'])[0].replace('_', ' ') if row['accepted_species'] else ''
                name = row['assigned_name'].replace('_', ' ') or 'unassigned'
                label = f"{name} | {row['readid']} | {folder}"
                metadata[identifier] = [identifier, 'local', name, '', row['readid'], folder, label]
                if not wanted:
                    break
        if wanted:
            raise ValueError(f'{folder}: tree read IDs missing from completed analysis: {list(wanted)[:3]}')
    return metadata


def render(tree_path, analysis_dir, reference_genbank, output_prefix):
    prefix = Path(output_prefix)
    outputs = [Path(str(prefix) + suffix) for suffix in
               ('.pdf', '.svg', '.midpoint.nwk', '.tips.tsv', '.plot.json')]
    for path in outputs:
        if path.exists():
            raise FileExistsError(path)
    tree_path = Path(tree_path)
    tree = bt.io.load_newick(str(tree_path), treeType='divergence')
    # Standard-bootstrap runs may put support only on the consensus tree.
    consensus = tree_path.with_suffix('.contree')
    if (tree_path.suffix == '.treefile' and consensus.is_file()
            and not any(b.traits.get('label') for b in tree.get_internal())):
        candidate = bt.io.load_newick(str(consensus), treeType='divergence')
        if any(b.traits.get('label') for b in candidate.get_internal()):
            if {b.name for b in candidate.get_external()} != {b.name for b in tree.get_external()}:
                raise ValueError('IQ-TREE consensus and ML trees have different tip sets')
            tree, tree_path = candidate, consensus
            print(f'Using supported consensus tree: {consensus}')
    ids = [leaf.name for leaf in tree.get_external()]
    if len(ids) != len(set(ids)) or len(ids) < 2:
        raise ValueError('Tree must contain at least two uniquely named tips')
    metadata = tip_metadata(ids, analysis_dir, reference_genbank)
    supports = midpoint_with_support(tree)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    # Vector figures remain usable when the tree has many tips.
    fig, ax = plt.subplots(figsize=(16, max(4, len(ids) * 0.22 + 1.5)))
    tree.plot_tree(ax, colour='black', width=0.8)
    span = max(leaf.x for leaf in tree.get_external())
    tree.plot_text(ax, targetFxn=lambda b: b.is_leaf(),
                   textContentFxn=lambda b: metadata[b.name][-1], xSpace=span * 0.012,
                   fontsize=8, va='center')
    for branch in tree.Objects:
        support = branch.traits.get('label')
        if support and branch is not tree.root:
            ax.annotate(support, ((branch.parent.x + branch.x) / 2, branch.y),
                        xytext=(0, 3), textcoords='offset points', ha='center', va='bottom', fontsize=7)
    bottom = min(leaf.y for leaf in tree.get_external())
    top = max(leaf.y for leaf in tree.get_external())
    plot_scale_bar(ax, xy=(0, bottom - max(1, tree.ySpan * 0.04)), tree=tree)
    ax.set_xlim(-span * 0.02, span * 1.05)
    ax.set_ylim(bottom - max(3, tree.ySpan * 0.12), top + 1)
    ax.set_axis_off()
    ax.set_title('Midpoint-rooted divergence tree', loc='left')
    for output in outputs[:2]:
        fig.savefig(output, bbox_inches='tight')
    plt.close(fig)
    outputs[2].write_text(tree.to_string(traits=['label']) + '\n')
    with outputs[3].open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['tree_tip_id', 'source', 'taxon_or_organism', 'accession', 'readid', 'source_folder', 'display_label'])
        writer.writerows(metadata[key] for key in ids)
    outputs[4].write_text(json.dumps({'baltic_version': version('baltic'), 'tree_type': 'divergence',
        'rooting': 'midpoint', 'tip_count': len(ids), 'supported_splits': len(supports),
        'support_labels': 'Unchanged IQ-TREE internal labels; slash-separated values retained',
        'scale_bar': 'Baltic default; substitutions/site', 'input_tree': str(Path(tree_path).resolve())}, indent=2) + '\n')
    print(f'Wrote {outputs[0]} and {outputs[1]} ({len(ids)} tips, {len(supports)} supported splits)')
    return tree


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tree', type=Path, required=True)
    parser.add_argument('--analysis-dir', type=Path, required=True)
    parser.add_argument('--reference-genbank', '--reference-metadata', dest='reference_genbank', type=Path, required=True, help='GenBank XML or BOLD reference metadata TSV')
    parser.add_argument('--output-prefix', type=Path, required=True)
    args = parser.parse_args()
    render(args.tree, args.analysis_dir, args.reference_genbank, args.output_prefix)


if __name__ == '__main__':
    main()
