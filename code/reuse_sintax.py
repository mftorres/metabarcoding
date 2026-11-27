#!/usr/bin/env python3
"""One-off for mergedseq5: classify the clusters, reusing the October 2026 read classifications.

mergedseq5 was classified read by read in October 2026 (full MIDORI database). After
re-cleaning and clustering it with the current Module 1 (CLEAN_ONLY=1), a cluster made
of a single read is that same read, so its October classification is copied. All
other clusters are classified now with SINTAX, with the database given by --db:
  * consensus sequences (clusters of two or more reads);
  * single reads that have no October row, or whose October pieces (fastplong splits
    chimeric reads) got different taxonomies.
The strand is flipped for reads that cutadapt reverse-complemented, so it refers to
the new read. Output: one classification file per sample, like a normal Module 1 run.
"""
import argparse
import csv
from pathlib import Path
import subprocess

from sample_folders import dataset_path, discover_samples, pipeline_setting
from seqio import fasta, raw_read_id, read_tsv, tsv, write_record

FLIP = {'+': '-', '-': '+'}


def old_classifications(path):
    """October SINTAX rows grouped by raw read ID: {raw ID: [(taxonomy, strand, cutoff column)]}."""
    rows = {}
    with open(path, encoding='utf-8', newline='') as handle:
        for fields in csv.reader(handle, delimiter='\t'):
            if fields:
                label = fields[3] if len(fields) > 3 else ''
                rows.setdefault(raw_read_id(fields[0].split()[0]), []).append((fields[1], fields[2], label))
    return rows


def classify(fasta_path, db, output):
    """Run SINTAX with the Module 1 settings."""
    subprocess.run(['vsearch', '--sintax', str(fasta_path), '--db', str(db), '--sintax_random', '--strand', 'both',
                    '--sintax_cutoff', pipeline_setting('SINTAX_CUTOFF'), '--threads', pipeline_setting('SINTAX_THREADS'),
                    '--randseed', pipeline_setting('SINTAX_RANDSEED'),
                    '--maxseqlength', pipeline_setting('SINTAX_MAX_LENGTH'),
                    '--minseqlength', pipeline_setting('SINTAX_MIN_LENGTH'), '--tabbedout', str(output)], check=True)


def reuse_folder(folder, old_path, db, output):
    filtered = folder/'filtered'
    old = old_classifications(old_path)
    flipped = {row['read_id'] for row in read_tsv(filtered/'primer_tags.tsv.gz') if row['reverse_complemented'] == 'yes'}
    sizes = {}
    for row in read_tsv(filtered/'cluster_members.tsv.gz'):
        sizes[row['cluster_id']] = sizes.get(row['cluster_id'], 0) + 1
    output.parent.mkdir(parents=True, exist_ok=True)
    reused = output.with_name(output.name + '.reused')
    to_classify = output.with_name(output.name + '.to_classify.fasta')
    classified = output.with_name(output.name + '.classified')
    counts = dict(consensus=0, single_reused=0, single_no_old_row=0, single_old_pieces_disagree=0)
    with reused.open('w') as copied, to_classify.open('w') as pending:
        for cluster_id, sequence in fasta(filtered/'clusters.fasta'):
            raw = raw_read_id(cluster_id)
            candidates = old.get(raw, [])
            if sizes[cluster_id] >= 2:
                counts['consensus'] += 1
            elif not candidates:
                counts['single_no_old_row'] += 1
            elif len({taxonomy for taxonomy, _, _ in candidates}) > 1:
                counts['single_old_pieces_disagree'] += 1
            else:
                counts['single_reused'] += 1
                taxonomy, strand, label = candidates[0]
                if raw in flipped:
                    strand = FLIP.get(strand, strand)
                copied.write(f'{cluster_id}\t{taxonomy}\t{strand}\t{label}\n')
                continue
            write_record(pending, cluster_id, sequence)
    classify(to_classify, db, classified)
    # Final file: reused rows followed by the new ones.
    with output.open('x') as handle:
        for part in (reused, classified):
            handle.write(part.read_text())
    for part in (reused, to_classify, classified):
        part.unlink()
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('dataset', help='Re-cleaned and clustered dataset, e.g. mergedseq5')
    parser.add_argument('--old-root', type=Path, required=True, help='Folder with the October <sample>/<sample>_SUFFIX.tsv files')
    parser.add_argument('--db', type=Path, required=True, help='SINTAX database (.udb) used in October')
    parser.add_argument('--suffix', default='sintaxSP', help='Suffix of the October files, also used for the new ones')
    parser.add_argument('--output-root', type=Path, help='Default data/sintax/DATASET')
    args = parser.parse_args()
    dataset = dataset_path(args.dataset)
    output_root = args.output_root or Path('data/sintax')/dataset.name
    if output_root.resolve() == args.old_root.resolve():
        parser.error('--output-root must differ from --old-root so the October results are kept')
    rows = []
    for folder in discover_samples([dataset]):
        old_path = args.old_root/folder.name/f'{folder.name}_{args.suffix}.tsv'
        output = output_root/folder.name/f'{folder.name}_{args.suffix}.tsv'
        if not old_path.is_file():
            print(f'[{folder.name}] no October classification; classify it with Module 1 instead', flush=True)
            continue
        if output.exists():
            print(f'[{folder.name}] {output} exists; skipped', flush=True)
            continue
        counts = reuse_folder(folder, old_path, args.db, output)
        rows.append([folder.name, *counts.values()])
        print(f'[{folder.name}] ' + ', '.join(f'{key} {value:,}' for key, value in counts.items()), flush=True)
    if rows:
        report = output_root/'reused_classification_report.tsv'
        if report.exists():
            rows = [list(row.values()) for row in read_tsv(report)] + rows
        tsv(report, ['source_folder', 'consensus', 'single_reused', 'single_no_old_row', 'single_old_pieces_disagree'], rows)


if __name__ == '__main__':
    main()
