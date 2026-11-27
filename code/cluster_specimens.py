#!/usr/bin/env python3
"""Cluster the cleaned reads of each specimen before taxonomic classification.

A specimen is one PCR well: one ONT barcode folder plus one forward tag. Reads
without an assigned forward tag are clustered together as specimen 'none'.
Within each specimen, VSEARCH groups reads that are at least 97% identical and
builds a consensus sequence for every cluster.

Outputs, next to the cleaned reads in FOLDER/filtered/:
  clusters.fasta           one sequence per cluster: the consensus when the cluster
                           has two or more reads, the read itself when it has one.
                           The sequence name (cluster ID) is the read ID of the
                           cluster's centroid (VSEARCH's seed: its longest read).
  cluster_members.tsv.gz   read_id, cluster_id and forward_tag of every cleaned read.
"""
import argparse
from pathlib import Path
import shutil
import subprocess
import tempfile

from seqio import fasta, fastq_records, raw_read_id, read_tsv, tsv, write_record


def forward_tags(tag_table):
    """Forward tag of every raw read; 'none' when the tag was not assigned."""
    tags = {}
    for row in read_tsv(tag_table):
        if row['forward_status'] == 'assigned':
            tags[row['read_id']] = row['forward_primer']
        else:
            tags[row['read_id']] = 'none'
    return tags


def split_by_specimen(reads, tags, work):
    """Write the reads of each specimen to its own FASTA file; return {specimen: path}."""
    paths, handles = {}, {}
    for read_id, sequence in fastq_records(reads):
        specimen = tags.get(raw_read_id(read_id))
        if specimen is None:
            raise ValueError(f'{read_id} has no row in the tag table; rerun the cleaning step')
        if specimen not in handles:
            paths[specimen] = work/f'{specimen}.fasta'
            handles[specimen] = paths[specimen].open('w')
        write_record(handles[specimen], read_id, sequence)
    for handle in handles.values():
        handle.close()
    return paths


def cluster(source, identity, threads, work):
    """Run VSEARCH on one specimen; return ({read: seed}, {seed: consensus sequence})."""
    stem = work/source.stem
    longest = max(len(sequence) for _, sequence in fasta(source))
    # Same settings as the reference clustering in Module 2: identity counted over the
    # aligned columns excluding terminal gaps (--iddef 2), both strands, no masking.
    subprocess.run(['vsearch', '--cluster_fast', str(source), '--id', str(identity), '--iddef', '2',
                    '--strand', 'both', '--qmask', 'none', '--minseqlength', '1',
                    '--maxseqlength', str(longest), '--notrunclabels', '--threads', str(threads),
                    '--consout', f'{stem}.consensus.fasta', '--uc', f'{stem}.uc', '--quiet'], check=True)
    seed_of = {}
    with open(f'{stem}.uc') as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            # S = seed of a new cluster, H = read added to the cluster of the seed in column 10.
            if fields[0] == 'S':
                seed_of[fields[8]] = fields[8]
            elif fields[0] == 'H':
                seed_of[fields[8]] = fields[9]
    consensus = {}
    for header, sequence in fasta(f'{stem}.consensus.fasta'):
        # VSEARCH names consensus sequences centroid=<seed>;seqs=<reads>.
        seed = header.removeprefix('centroid=').rpartition(';seqs=')[0]
        consensus[seed] = sequence
    return seed_of, consensus


def cluster_folder(folder, identity, threads):
    filtered = Path(folder)/'filtered'
    tags = forward_tags(filtered/'primer_tags.tsv.gz')
    work = Path(tempfile.mkdtemp(prefix='.clusters.', dir=filtered))
    try:
        paths = split_by_specimen(filtered/'asinput.fastq.gz', tags, work)
        members, consensus_count, singleton_count = [], 0, 0
        with (work/'clusters.fasta').open('w') as output:
            for specimen in sorted(paths):
                seed_of, consensus = cluster(paths[specimen], identity, threads, work)
                sizes = {}
                for read_id, seed in seed_of.items():
                    sizes[seed] = sizes.get(seed, 0) + 1
                    members.append((read_id, seed, specimen))
                reads = dict(fasta(paths[specimen]))
                for seed in sorted(sizes):
                    # A one-read cluster keeps the read itself, not VSEARCH's copy of it.
                    if sizes[seed] >= 2:
                        write_record(output, seed, consensus[seed])
                        consensus_count += 1
                    else:
                        write_record(output, seed, reads[seed])
                        singleton_count += 1
        tsv(work/'cluster_members.tsv.gz', ['read_id', 'cluster_id', 'forward_tag'], members)
        # Publish both files only when everything succeeded.
        for name in ('cluster_members.tsv.gz', 'clusters.fasta'):
            (work/name).replace(filtered/name)
    finally:
        shutil.rmtree(work)
    print(f'[{Path(folder).name}] {len(members):,} cleaned reads in {len(paths)} specimens -> '
          f'{consensus_count:,} consensus sequences + {singleton_count:,} single reads '
          f'= {consensus_count + singleton_count:,} sequences to classify', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('folder', help='Sample folder (barcodeNN or unclassified) with a filtered/ subfolder')
    parser.add_argument('--identity', type=float, default=0.97)
    parser.add_argument('--threads', type=int, default=8)
    args = parser.parse_args()
    cluster_folder(args.folder, args.identity, args.threads)


if __name__ == '__main__':
    main()
