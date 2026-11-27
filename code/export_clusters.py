#!/usr/bin/env python3
"""Export cluster sequences (consensus sequences and single reads) from a Module 1 analysis.

Select clusters by taxon (--taxon: accepted with confidence >= the analysis cutoff at
any rank), by forward tag (--tag), or both. Sequences come from each sample's
filtered/clusters.fasta. Each exported sequence is named

    folder-forwardtag-taxon-contaminationflag-readcount-clusterID

Fields are percent-encoded so they never contain '-' (the cluster ID, last, may).
"""
import argparse
import difflib
import json
from pathlib import Path
import re
from urllib.parse import quote, unquote

from sample_folders import analysis_folders, from_project
from seqio import fasta, read_tsv, write_record
from sintax_parser import RANKS, taxon_parts

FLAGS = ('contamination', 'nocontamination')


def header_field(value):
    return quote(re.sub(r'\s+', '_', value), safe='_.').replace('-', '%2D')


def export_header(folder, forward_tag, taxon, contamination, read_count, cluster_id):
    fields = [folder, forward_tag, taxon or 'unassigned']
    return '-'.join([header_field(value) for value in fields] + [contamination, str(read_count), cluster_id])


def header_fields(identifier):
    """(folder, forward tag, taxon, flag, read count, cluster ID) of an exported cluster, or None
    for other sequences (references)."""
    parts = identifier.split()[0].split('-', 5)
    if len(parts) != 6 or parts[3] not in FLAGS or not parts[4].isdigit():
        return None
    return unquote(parts[0]), unquote(parts[1]), unquote(parts[2]), parts[3], int(parts[4]), parts[5]


def normalise(value):
    return re.sub(r'[\s_]+', ' ', value.strip()).casefold()


def analysis_inputs(analysis_dir, folder):
    """Cluster table and cluster sequences of one folder; refuse sequences changed since the analysis."""
    analysis_dir = Path(analysis_dir)
    if not (analysis_dir/'analysis.complete').is_file():
        raise ValueError(f'Analysis is incomplete: {analysis_dir}')
    params = json.loads((analysis_dir/'analysis_parameters.json').read_text())
    if 'cluster_identity' not in params:
        raise ValueError(f'{analysis_dir} was made before clustering moved to Module 1 (November 2026); '
                         'rerun Module 1 on the dataset and use the new analysis folder')
    record = params['input_files'][folder]['clusters']
    clusters = from_project(record['path'])
    if clusters.stat().st_size != record['bytes']:
        raise ValueError(f'{clusters} changed after the analysis; rerun the Module 1 analysis first')
    return analysis_dir/folder/f'{folder}_clusters.tsv', clusters


def matches(row, targets, tags):
    """True when the cluster has one of the forward tags and is accepted as one of the taxa."""
    if tags and row['forward_tag'] not in tags:
        return False
    if not targets:
        return True
    for rank in RANKS.values():
        if row[rank]:
            name, taxid = taxon_parts(row[rank])
            if targets & {normalise(name), normalise(row[rank]), taxid}:
                return True
    return False


def export_clusters(*, analysis_dir, output, taxa=(), tags=(), folders=None):
    """Write the selected clusters to output; return the number written."""
    output = Path(output)
    if output.exists():
        raise FileExistsError(f'Refusing to overwrite {output}')
    if not taxa and not tags:
        raise ValueError('Give at least one --taxon or --tag')
    targets = {normalise(taxon) for taxon in taxa}
    output.parent.mkdir(parents=True, exist_ok=True)
    # Write to a .partial file and give it its final name only when every folder succeeded.
    partial = output.with_name(output.name + '.partial')
    total, seen_names = 0, set()
    try:
        with partial.open('w') as handle:
            for folder in analysis_folders(folders, analysis_dir):
                table, sequences = analysis_inputs(analysis_dir, folder)
                selected = {}
                for row in read_tsv(table):
                    if row['assigned_name']:
                        seen_names.add(normalise(row['assigned_name']))
                    if matches(row, targets, set(tags)):
                        selected[row['cluster_id']] = export_header(folder, row['forward_tag'], row['assigned_label'],
                                                                    row['contamination'], row['read_count'],
                                                                    row['cluster_id'])
                count = 0
                for cluster_id, sequence in fasta(sequences):
                    if cluster_id in selected:
                        write_record(handle, selected[cluster_id], sequence)
                        count += 1
                if count != len(selected):
                    raise ValueError(f'{folder}: {len(selected) - count} selected clusters are missing from {sequences}')
                total += count
                print(f'{folder}: exported {count:,} clusters', flush=True)
        if not total:
            hint = {taxon: difflib.get_close_matches(normalise(taxon), seen_names, n=4) for taxon in taxa}
            raise ValueError(f'No clusters match taxa {list(taxa)} and tags {list(tags)}. Closest taxon names: {hint}')
    except BaseException:
        partial.unlink(missing_ok=True)
        raise
    partial.replace(output)
    print(f'Wrote {total:,} sequences to {output}', flush=True)
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--analysis-dir', type=Path, required=True, help='Module 1 analysis folder')
    parser.add_argument('--output', type=Path, required=True, help='New FASTA file')
    parser.add_argument('--taxon', dest='taxa', action='append', default=[],
                        help='Taxon name, MIDORI label or NCBI taxid; repeat for several')
    parser.add_argument('--tag', dest='tags', action='append', default=[], help='Forward tag, e.g. bc1087; repeat for several')
    parser.add_argument('--folders', nargs='+', help='Sample folders; default every analysed folder')
    args = parser.parse_args()
    try:
        export_clusters(analysis_dir=args.analysis_dir, output=args.output, taxa=args.taxa, tags=args.tags,
                        folders=args.folders)
    except (OSError, ValueError) as error:
        parser.exit(1, f'Export failed: {error}\n')


if __name__ == '__main__':
    main()
