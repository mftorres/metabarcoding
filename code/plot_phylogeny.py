#!/usr/bin/env python3
"""Render an IQ-TREE divergence phylogeny with Baltic 1.x and midpoint rooting."""
import argparse
import csv
from importlib.metadata import version
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import baltic as bt
from baltic.bt_utils import plot_scale_bar
from export_clusters import header_fields
from reference_labels import reference_names, reference_id
from sintax_parser import taxon_parts


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


def tip_metadata(tip_ids, reference_genbank):
    """Describe each tip: local clusters from their sequence names, references from their metadata."""
    references = reference_names(reference_genbank)
    bold_codes = {}
    if Path(reference_genbank).suffix == '.tsv':
        with Path(reference_genbank).open() as handle:
            bold_codes = {row['sequence_id']: row['process_id'] for row in csv.DictReader(handle, delimiter='\t')}
    metadata = {}
    for identifier in tip_ids:
        fields = header_fields(identifier)
        if fields is not None:
            folder, forward, taxon, _, read_count, cluster_id = fields
            name = taxon_parts(taxon)[0].replace('_', ' ')
            label = f'{name} | {folder} | {forward} | {read_count} reads'
            metadata[identifier] = [identifier, 'local', name, '', cluster_id, read_count, folder, forward, label]
        else:
            # efetch FASTA IDs may be accession.version or pipe-separated.
            base_id = reference_id(identifier)
            candidates = [base_id] + base_id.split('|')
            accession = next((x for x in candidates if x in references), None)
            if accession is None:
                raise ValueError(f'No source metadata for reference {identifier}')
            organism = references[accession]
            display_id = bold_codes.get(accession, accession)
            metadata[identifier] = [identifier, 'bold' if bold_codes else 'reference', organism, display_id,
                                    '', '', '', '', f'{organism} | {display_id}']
    return metadata


def render(tree_path, reference_genbank, output_prefix):
    prefix = Path(output_prefix)
    outputs = [Path(str(prefix) + suffix) for suffix in
               ('.pdf', '.svg', '.midpoint.nwk', '.tips.tsv', '.plot.json', '.png')]
    for path in outputs:
        if path.exists():
            raise FileExistsError(path)
    tree = bt.io.load_newick(str(tree_path), treeType='divergence')
    ids = [leaf.name for leaf in tree.get_external()]
    if len(ids) != len(set(ids)) or len(ids) < 2:
        raise ValueError('Tree must contain at least two uniquely named tips')
    metadata = tip_metadata(ids, reference_genbank)
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
    for output in [*outputs[:2], outputs[5]]:
        fig.savefig(output, bbox_inches='tight', dpi=100)
    plt.close(fig)
    outputs[2].write_text(tree.to_string(traits=['label']) + '\n')
    with outputs[3].open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['tree_tip_id', 'source', 'taxon_or_organism', 'accession', 'cluster_id', 'read_count',
                         'source_folder', 'forward_tag', 'display_label'])
        writer.writerows(metadata[key] for key in ids)
    outputs[4].write_text(json.dumps({'baltic_version': version('baltic'), 'tree_type': 'divergence',
        'rooting': 'midpoint', 'tip_count': len(ids), 'supported_splits': len(supports),
        'support_labels': 'Unchanged IQ-TREE internal labels; slash-separated values retained',
        'scale_bar': 'Baltic default; substitutions/site', 'input_tree': str(Path(tree_path).resolve())}, indent=2) + '\n')
    print(f'Wrote {outputs[0]}, {outputs[1]} and {outputs[5]} ({len(ids)} tips, {len(supports)} supported splits)')
    return tree


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tree', type=Path, required=True)
    parser.add_argument('--reference-genbank', '--reference-metadata', dest='reference_genbank', type=Path, required=True, help='GenBank XML or BOLD reference metadata TSV')
    parser.add_argument('--output-prefix', type=Path, required=True)
    args = parser.parse_args()
    render(args.tree, args.reference_genbank, args.output_prefix)


if __name__ == '__main__':
    main()
