#!/usr/bin/env python3
"""Export annotated cleaned reads whose raw parents have a specified best forward tag."""
import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

from export_sintax_fasta import fastq_records
from sequence_headers import export_header
from sintax_parser import parse_output


def clean_records(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence, plus, quality = handle.readline().rstrip('\r\n'), handle.readline(), handle.readline().rstrip('\r\n')
            if not header.startswith('@') or not plus.startswith('+') or len(sequence) != len(quality):
                raise ValueError(f'Invalid FASTQ in {path}')
            readid = header[1:].split()[0]
            parent = re.search(r'(?:^|\s)parent_read_id=([^\s]+)', header)
            raw_id = parent[1] if parent else re.sub(r'^split-by-adapter-(?:left|right)-', '', readid)
            yield readid, raw_id


def select_hits(path, dataset, folders, primer, mapping):
    canonical = {value for values in mapping.values() for value in values.split(';')}
    if primer not in canonical:
        raise ValueError(f'Unknown canonical primer {primer!r}; use a name from primer_name_mapping.json')
    selected = {folder: {} for folder in folders}
    seen = {folder: set() for folder in folders}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'merged_folder', 'barcode_folder', 'seqID', 'forward_patternName'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f'Best-hit table is missing columns: {required - set(reader.fieldnames or [])}')
        for row in reader:
            folder = row['barcode_folder']
            if row['merged_folder'] != dataset or folder not in selected:
                continue
            raw_id = row['seqID']
            if raw_id in seen[folder]:
                raise ValueError(f'Duplicate best-hit row for {folder}/{raw_id}')
            seen[folder].add(raw_id)
            original = row['forward_patternName']
            if not original:
                continue
            mapped = mapping.get(original, original)
            if not set(mapped.split(';')).issubset(canonical):
                raise ValueError(f'Unmapped forward-primer name: {original}')
            if primer in mapped.split(';'):
                reverse = row.get('reverse_patternName', '')
                reverse = mapping.get(reverse, reverse)
                if reverse and not set(reverse.split(';')).issubset(canonical):
                    raise ValueError(f'Unmapped reverse-primer name: {reverse}')
                selected[folder][raw_id] = (mapped, reverse or 'none')
    for folder in folders:
        if not seen[folder]:
            raise ValueError(f'No best-hit tracking rows for {dataset}/{folder}')
    return selected


def raw_files(folder):
    paths = sorted({p for pattern in ('*.fastq', '*.fastq.gz', '*.fq', '*.fq.gz') for p in folder.glob(pattern)})
    if not paths:
        raise FileNotFoundError(f'No top-level raw FASTQ files in {folder}')
    return paths


def export_primer(*, primer, barcode, analysis_dir, output, include_unclassified=False,
                  primer_best_hits=None, merged_dir=None, threads=4):
    if not shutil.which('seqkit'):
        raise ValueError('seqkit is not on PATH; activate the taxonomy environment')
    analysis_dir, output = Path(analysis_dir), Path(output)
    if output.suffix.lower() not in ('.fasta', '.fa', '.fna'):
        raise ValueError('--output must be an uncompressed .fasta, .fa or .fna file')
    if not (analysis_dir / 'analysis.complete').is_file():
        raise ValueError('Provide a completed SINTAX analysis directory')
    params = json.loads((analysis_dir / 'analysis_parameters.json').read_text())
    if params['cutoff'] != 0.8:
        raise ValueError('Headers require a completed analysis with cutoff 0.80')
    mapping = json.loads((analysis_dir / 'primer_name_mapping.json').read_text())
    merged_dir = Path(merged_dir or params['merged_dir'])
    primer_best_hits = Path(primer_best_hits or params['primer_hits'])
    folders = [Path(barcode).name]
    if include_unclassified and 'unclassified' not in folders:
        folders.append('unclassified')
    for folder in folders:
        if folder not in params['analysed_folders']:
            raise ValueError(f'{folder} is not in this completed analysis. Clean/classify it and include it in the analysis before exporting annotated reads.')
    destinations = [output, output.with_suffix('.reads.tsv'), output.with_suffix('.counts.tsv'), output.with_suffix('.log')]
    for path in destinations:
        if path.exists():
            raise FileExistsError(f'Refusing to overwrite {path}')
    with (analysis_dir / 'input_manifest.tsv').open() as handle:
        manifest = {row['source_folder']: row for row in csv.DictReader(handle, delimiter='\t')}
    with (analysis_dir / 'taxon_plate_prevalence.tsv').open() as handle:
        flagged = {(row['rank'], row['taxon_label']) for row in csv.DictReader(handle, delimiter='\t')
                   if row['contamination'] == 'contamination'}
    sources = {}
    for folder in folders:
        source = Path(manifest[folder]['sintax_file'])
        stat = source.stat()
        if str(stat.st_size) != manifest[folder]['sintax_bytes'] or str(stat.st_mtime_ns) != manifest[folder]['sintax_mtime_ns']:
            raise ValueError(f'SINTAX changed since analysis: {source}')
        clean = merged_dir / folder / 'filtered/asinput.fastq.gz'
        if not clean.is_file():
            raise FileNotFoundError(clean)
        sources[folder] = (source, clean, raw_files(merged_dir / folder))
    selected = select_hits(primer_best_hits, merged_dir.name, folders, primer, mapping)
    if not any(selected.values()):
        raise ValueError(f'No raw reads have best forward tag {primer} in {", ".join(folders)}')
    output.parent.mkdir(parents=True, exist_ok=True)
    total = 0
    with tempfile.TemporaryDirectory(prefix='.primer_export_', dir=output.parent) as directory:
        temp = Path(directory)
        staged = [temp / p.name for p in destinations]
        with staged[0].open('w') as fasta_out, staged[1].open('w', newline='') as audit, staged[3].open('w') as log:
            fields = ['source_folder', 'forward_primers_detected', 'reverse_primers_detected', 'raw_ont_read_id', 'raw_fastq_file',
                      'readid', 'clean_fastq_file', 'assigned_rank', 'assigned_label', 'assigned_confidence',
                      'taxonomy_status', 'contamination', 'export_header', 'status']
            writer = csv.DictWriter(audit, fields, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            log.write(f'primer={primer}\nbest_hits={primer_best_hits.resolve()}\nanalysis={analysis_dir.resolve()}\n')
            counts = []
            def grep_ids(ids, source, target):
                if not ids:
                    raise ValueError('Internal error: refusing SeqKit grep with an empty ID list')
                patterns = temp / 'ids.txt'
                patterns.write_text(''.join(identifier + '\n' for identifier in sorted(ids)))
                command = ['seqkit', 'grep', '--threads', str(threads), '--pattern-file', str(patterns), str(source)]
                log.write(' '.join(command) + '\n'); log.flush()
                with target.open('w') as handle:
                    subprocess.run(command, stdout=handle, stderr=log, check=True)
            for folder in folders:
                hits = selected[folder]
                source, clean, paths = sources[folder]
                if not hits:
                    counts.append((folder, 0, 0, 0, 0))
                    continue
                parents = {}
                for raw_path in paths:
                    matches = temp / 'raw_matches.fastq'
                    grep_ids(hits, raw_path, matches)
                    for raw_id, _ in fastq_records(matches):
                        if raw_id in parents:
                            raise ValueError(f'Duplicate raw read ID in {folder}: {raw_id}')
                        parents[raw_id] = str(raw_path.resolve())
                if set(parents) != set(hits):
                    raise ValueError(f'{folder}: {len(set(hits)-set(parents))} primer-selected IDs missing from raw FASTQ files')
                wanted = {}
                for readid, raw_id in clean_records(clean):
                    if raw_id in hits:
                        if readid in wanted:
                            raise ValueError(f'Duplicate cleaned read ID: {readid}')
                        wanted[readid] = raw_id
                assignments = {}
                for row in parse_output(source, cutoff=params['cutoff']):
                    if row['readid'] in wanted:
                        if row['readid'] in assignments or row['raw_ont_read_id'] != wanted[row['readid']]:
                            raise ValueError(f'Duplicate or inconsistent SINTAX parent mapping: {row["readid"]}')
                        assignments[row['readid']] = row
                if set(assignments) != set(wanted):
                    raise ValueError(f'{folder}: {len(set(wanted)-set(assignments))} selected cleaned reads lack SINTAX rows; complete classification/analysis first')
                matched_parents = set(wanted.values())
                if wanted:
                    matches = temp / 'clean_matches.fastq'
                    grep_ids(wanted, clean, matches)
                    for readid, sequence in fastq_records(matches):
                        raw_id = wanted.pop(readid)
                        row = assignments[readid]
                        flag = 'contamination' if (row['assigned_rank'], row['assigned_label']) in flagged else 'nocontamination'
                        header = export_header(folder, *hits[raw_id], row['assigned_label'], flag, readid)
                        fasta_out.write('>' + header + '\n')
                        for offset in range(0, len(sequence), 80):
                            fasta_out.write(sequence[offset:offset+80] + '\n')
                        writer.writerow(dict(source_folder=folder, forward_primers_detected=hits[raw_id][0], reverse_primers_detected=hits[raw_id][1],
                            raw_ont_read_id=raw_id, raw_fastq_file=parents[raw_id], readid=readid,
                            clean_fastq_file=str(clean.resolve()), assigned_rank=row['assigned_rank'],
                            assigned_label=row['assigned_label'], assigned_confidence=row['assigned_confidence'],
                            taxonomy_status=row['taxonomy_status'], contamination=flag, export_header=header, status='exported'))
                        total += 1
                    if wanted:
                        raise ValueError(f'{folder}: SeqKit did not return all selected cleaned reads')
                for raw_id in sorted(set(hits)-matched_parents):
                    writer.writerow(dict(source_folder=folder, forward_primers_detected=hits[raw_id][0], reverse_primers_detected=hits[raw_id][1],
                        raw_ont_read_id=raw_id, raw_fastq_file=parents[raw_id], clean_fastq_file=str(clean.resolve()),
                        status='no_clean_read'))
                counts.append((folder, len(hits), len(matched_parents), len(assignments), len(hits)-len(matched_parents)))
                message = f'{folder}: {len(hits)} raw best-tag hits; {len(assignments)} cleaned reads exported; {len(hits)-len(matched_parents)} raw parents have no clean reads'
                print(message, flush=True); log.write(message + '\n')
            log.write(f'Total cleaned reads exported: {total}\n')
        with staged[2].open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(['source_folder', 'raw_best_tag_hits', 'raw_parents_with_clean_reads', 'clean_reads_exported', 'raw_parents_without_clean_reads'])
            writer.writerows(counts)
            writer.writerow(['TOTAL'] + [sum(row[i] for row in counts) for i in range(1,5)])
        # Publish exclusively after all folders validate; never replace existing outputs.
        published = []
        try:
            for stage, destination in zip(staged, destinations):
                os.link(stage, destination)
                published.append(destination)
        except BaseException:
            for path in published:
                path.unlink()
            raise
    print(f'Wrote {total} cleaned reads to {output}')
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--primer', required=True, help='Canonical forward tag, e.g. bc1087')
    parser.add_argument('--barcode', required=True, help='Folder name, e.g. barcode15')
    parser.add_argument('--analysis-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--include-unclassified', action='store_true', help='Also search the unclassified folder (off by default)')
    parser.add_argument('--primer-best-hits', type=Path, help='Defaults to the best-hit table recorded in the analysis')
    parser.add_argument('--merged-dir', type=Path, help='Defaults to the raw/clean dataset recorded in the analysis')
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error('--threads must be positive')
    try:
        export_primer(**vars(args))
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        parser.exit(1, f'Primer export failed: {error}\n')


if __name__ == '__main__':
    main()
