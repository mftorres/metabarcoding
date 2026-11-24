#!/usr/bin/env python3
"""Copy classifications from an earlier SINTAX run onto newly cleaned reads, matched by raw read ID.

Use this when a dataset was re-cleaned with the current Module 1 (CLEAN_ONLY=1) but
its reads were already classified by an older pipeline version, so the weeks-long
SINTAX step need not be repeated. The old classifications were made on reads that
still carried primers and tags; the classified COI region is the same.

Matching, per cleaned read:
  * the raw read has exactly one old SINTAX row -> copied;
  * the raw read was split into several old pieces that all got the same
    taxonomy -> copied;
  * pieces with different taxonomies, or no old row -> left out (counted in the report).
The strand is flipped for reads that cutadapt reverse-complemented, so it stays
relative to the new read.
"""
import argparse
import csv
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path

from sample_folders import dataset_path, discover_samples, project_path
from seqio import fastq_records
from sintax_parser import raw_read_id

FLIP = {'+': '-', '-': '+'}
REPORT_COLUMNS = ['source_folder', 'cleaned_reads', 'copied_single', 'copied_split_agree',
                  'skipped_split_disagree', 'skipped_not_in_old', 'old_rows', 'old_raw_reads_unused']


def old_rows(path):
    rows = {}
    with path.open(encoding='utf-8', newline='') as handle:
        for fields in csv.reader(handle, delimiter='\t'):
            if fields:
                header, taxonomy, strand = fields[:3]
                label = fields[3] if len(fields) > 3 else ''
                rows.setdefault(raw_read_id(header.split()[0]), []).append((taxonomy, strand, label))
    return rows


def reversed_reads(path):
    with gzip.open(path, 'rt', newline='') as handle:
        return {row['read_id'] for row in csv.DictReader(handle, delimiter='\t')
                if row['reverse_complemented'] == 'yes'}


def reuse_folder(folder, old_path, output):
    old = old_rows(old_path)
    flipped = reversed_reads(folder/'filtered/primer_tags.tsv.gz')
    counts = dict.fromkeys(REPORT_COLUMNS[1:], 0)
    counts['old_rows'] = sum(map(len, old.values()))
    used = set()
    temporary = output.with_suffix('.tmp')
    output.parent.mkdir(parents=True, exist_ok=True)
    with temporary.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        for readid, _ in fastq_records(folder/'filtered/asinput.fastq.gz'):
            counts['cleaned_reads'] += 1
            raw = raw_read_id(readid)
            candidates = old.get(raw)
            if not candidates:
                counts['skipped_not_in_old'] += 1
                continue
            if len({c[0] for c in candidates}) > 1:
                counts['skipped_split_disagree'] += 1
                continue
            counts['copied_single' if len(candidates) == 1 else 'copied_split_agree'] += 1
            used.add(raw)
            taxonomy, strand, label = candidates[0]
            writer.writerow([readid, taxonomy, FLIP.get(strand, strand) if raw in flipped else strand, label])
    temporary.replace(output)
    counts['old_raw_reads_unused'] = len(old) - len(used)
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('dataset', help='Re-cleaned dataset, e.g. mergedseq5')
    parser.add_argument('--old-root', type=Path, required=True,
                        help='Folder holding the old <sample>/<sample>_SUFFIX.tsv files, e.g. data/sintax')
    parser.add_argument('--output-root', type=Path, help='New classification folder; default data/sintax/DATASET')
    parser.add_argument('--suffix', default='sintaxSP')
    args = parser.parse_args()
    dataset = dataset_path(args.dataset)
    output_root = args.output_root or Path('data/sintax')/dataset.name
    if output_root.resolve() == args.old_root.resolve():
        parser.error('--output-root must differ from --old-root so the old results are kept')
    rows, sources = [], {}
    for folder in discover_samples([dataset]):
        old_path = args.old_root/folder.name/f'{folder.name}_{args.suffix}.tsv'
        if not old_path.is_file():
            print(f'[{folder.name}] no old classification; skipped')
            continue
        output = output_root/folder.name/f'{folder.name}_{args.suffix}.tsv'
        if output.exists():
            parser.error(f'{output} exists; remove it first')
        counts = reuse_folder(folder, old_path, output)
        rows.append([folder.name, *counts.values()])
        sources[folder.name] = {'old_sintax': project_path(old_path), 'old_bytes': old_path.stat().st_size,
                                'output': project_path(output)}
        copied = counts['copied_single'] + counts['copied_split_agree']
        print(f'[{folder.name}] copied {copied}/{counts["cleaned_reads"]} '
              f'({100*copied/max(counts["cleaned_reads"], 1):.1f}%)', flush=True)
    if not rows:
        parser.error('No sample had an old classification')
    with (output_root/'reused_classification_report.tsv').open('x', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerows([REPORT_COLUMNS, *rows])
    (output_root/'reused_classification.json').write_text(json.dumps({
        'created_utc': datetime.now(timezone.utc).isoformat(), 'method': 'old SINTAX rows matched by raw read ID',
        'note': 'Old classifications were made on reads that still carried primers and tags.',
        'samples': sources}, indent=2) + '\n')


if __name__ == '__main__':
    main()
