#!/usr/bin/env python3
"""Export cleaned reads confidently assigned (>=0.80 at the queried rank) to a taxon.

Selection and annotation come from the completed analysis read tables; sequences come
from each sample's cleaned, primer-trimmed reads (filtered/asinput.fastq.gz).
"""
import argparse
import csv
import difflib
import json
import os
from pathlib import Path
import re
import tempfile

from sample_folders import analysis_folders, from_project
from sequence_headers import export_header
from seqio import fastq_records, write_record
from sintax_parser import RANKS, taxon_parts


def normalize(value):
    return re.sub(r'[\s_]+', ' ', value.strip()).casefold()


def analysis_inputs(analysis_dir, folder):
    """Return a folder's read table and cleaned reads, refusing files changed since the analysis."""
    analysis_dir = Path(analysis_dir)
    if not (analysis_dir/'analysis.complete').is_file():
        raise ValueError(f'Analysis is incomplete: {analysis_dir}')
    params = json.loads((analysis_dir/'analysis_parameters.json').read_text())
    if folder not in params['analysed_folders']:
        raise ValueError(f'{folder} was not analysed in {analysis_dir}')
    if 'input_files' not in params:
        raise ValueError(f'{analysis_dir} was made by an older pipeline version (before Nov 2026) and does not '
                         'record its cleaned reads; rerun sintax_analysis.py on the dataset and use the new analysis folder')
    record = params['input_files'][folder]['clean_reads']
    clean = from_project(record['path'])
    if clean.stat().st_size != record['bytes']:
        raise ValueError(f'{clean} changed after the analysis; rerun the Read processing analysis first')
    return analysis_dir/folder/f'{folder}_reads.tsv', clean


def read_rows(table):
    with Path(table).open(newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def write_selected(writer, clean, selected, folder):
    """Copy selected reads (readid -> header) from the cleaned FASTQ; return the count."""
    wanted, count = dict(selected), 0
    for readid, sequence in fastq_records(clean):
        header = wanted.pop(readid, None)
        if header is not None:
            write_record(writer, header, sequence)
            count += 1
    if wanted:
        raise ValueError(f'{folder}: {len(wanted)} selected reads absent from {clean}, e.g. {list(wanted)[:3]}')
    return count


def export_reads(*, taxa, folders, analysis_dir, output):
    output = Path(output)
    if output.exists():
        raise FileExistsError(output)
    targets = {normalize(t) for t in taxa}
    seen_taxa, total = set(), 0
    output.parent.mkdir(parents=True, exist_ok=True)
    # Build privately and publish only after every folder has been exported.
    with tempfile.TemporaryDirectory(prefix='.sintax_fasta_', dir=output.parent) as temp:
        temporary = Path(temp)/'reads.fasta'
        with temporary.open('w') as writer:
            for folder in analysis_folders(folders, analysis_dir):
                table, clean = analysis_inputs(analysis_dir, folder)
                selected = {}
                for row in read_rows(table):
                    labels = [row['accepted_' + rank] for rank in RANKS.values() if row['accepted_' + rank]]
                    keys = set()
                    for label in labels:
                        name, taxid = taxon_parts(label)
                        seen_taxa.add(normalize(name))
                        keys |= {normalize(name), normalize(label), taxid}
                    if targets & keys:
                        selected[row['readid']] = export_header(folder, row['forward_PCR_primer'], row['reverse_PCR_primer'],
                                                                row['assigned_label'], row['contamination'], row['readid'])
                count = write_selected(writer, clean, selected, folder)
                total += count
                print(f'{folder}: exported {count:,} cleaned reads', flush=True)
        if not total:
            suggestions = {t: difflib.get_close_matches(normalize(t), seen_taxa, n=4) for t in taxa}
            raise ValueError(f'No reads assigned to {taxa} with confidence >=0.80. Closest names: {suggestions}')
        os.link(temporary, output)  # exclusive publication; never replace an existing FASTA
    print(f'Wrote {total:,} sequences to {output}', flush=True)
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--taxon', dest='taxa', action='append', required=True,
                        help='Taxon name, MIDORI label or NCBI taxid; repeat to export a union')
    parser.add_argument('--folders', nargs='+', help='Sample folders; default: every analysed folder')
    parser.add_argument('--analysis-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        export_reads(taxa=args.taxa, folders=args.folders, analysis_dir=args.analysis_dir, output=args.output)
    except (OSError, ValueError) as error:
        parser.exit(1, f'Export failed: {error}\n')


if __name__ == '__main__':
    main()
