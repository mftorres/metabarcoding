#!/usr/bin/env python3
"""Read-level SINTAX analysis, primer tags and barcode prevalence for one dataset."""
import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import tempfile

import pandas as pd
from sintax_parser import RANKS, COLUMNS, parse_output, taxon_parts
from sample_folders import discover_samples, dataset_path, project_path
from analysis_report import write_report
from describe_table_headers import describe_headers

RANK_NAMES = list(RANKS.values())
TAG_FIELDS = {'forward_primer': 'forward_PCR_primer', 'reverse_primer': 'reverse_PCR_primer',
              'forward_status': 'forward_tag_status', 'reverse_status': 'reverse_tag_status'}


def join_values(values):
    return ';'.join(sorted({name for v in values if pd.notna(v) for name in str(v).split(';') if name}))


def write_tsv(frame, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep='\t', index=False, header=True, mode='x', na_rep='')


def file_record(path):
    return dict(path=project_path(path), bytes=Path(path).stat().st_size)


def sample_frame(path, cutoff, tags):
    """Parse one SINTAX file and attach the tag calls of each read's raw parent."""
    frame = pd.DataFrame.from_records(parse_output(path, cutoff=cutoff), columns=COLUMNS)
    if frame.readid.duplicated().any():
        raise ValueError(f'Duplicate query IDs in {path}')
    for col in ['assigned_confidence', 'taxonomy_lowestboot'] + [r + '_bootstrap' for r in RANK_NAMES]:
        frame[col] = pd.to_numeric(frame[col], errors='coerce')
    missing = ~frame.raw_ont_read_id.isin(tags.index)
    if missing.any():
        raise ValueError(f'{path}: {missing.sum()} reads have no primer_tags.tsv.gz row '
                         f'(e.g. {frame.raw_ont_read_id[missing].iloc[0]}); rerun Read processing for this sample')
    for source, target in TAG_FIELDS.items():
        frame[target] = frame.raw_ont_read_id.map(tags[source]).fillna('')
    return frame


def prevalence_table(presence, barcodes, fraction):
    minimum = math.ceil(fraction * len(barcodes) - 1e-12)
    rows = []
    for (rank, label), observed in sorted(presence.items()):
        name, taxid = taxon_parts(label)
        rows.append(dict(rank=rank, taxon_label=label, taxon_name=name, taxid=taxid,
                         barcode_count=len(observed), analysed_barcode_count=len(barcodes),
                         barcode_fraction=len(observed)/len(barcodes), minimum_barcode_count=minimum,
                         prevalence_threshold=fraction, barcodes=';'.join(sorted(observed)),
                         contamination='contamination' if len(observed) >= minimum else 'nocontamination'))
    return pd.DataFrame(rows, columns=['rank', 'taxon_label', 'taxon_name', 'taxid', 'barcode_count',
        'analysed_barcode_count', 'barcode_fraction', 'minimum_barcode_count', 'prevalence_threshold',
        'barcodes', 'contamination'])


def add_flags(frame, flagged):
    frame['contamination'] = [
        'contamination' if (rank, label) in flagged else 'nocontamination'
        for rank, label in zip(frame.assigned_rank, frame.assigned_label)]
    for rank in RANK_NAMES:
        frame[rank + '_contamination'] = frame['accepted_' + rank].map(
            lambda label: 'contamination' if (rank, label) in flagged else 'nocontamination')
    return frame


def species_summary(frame, flagged, cutoff):
    """One row per candidate species; statistics distinguish confident and candidate calls."""
    candidates = frame[frame.species.ne('')].copy()
    keys = ['source_folder', 'species']
    columns = keys + ['species_name', 'species_taxid', 'genus', 'family', 'order', 'class',
                      'contamination', 'median_confidence', 'read_count', 'raw_read_count',
                      'forward_primers_detected', 'reverse_primers_detected',
                      'candidate_read_count', 'median_candidate_confidence']
    if candidates.empty:
        return pd.DataFrame(columns=columns)
    all_species = candidates.groupby(keys, dropna=False).agg(
        genus=('genus', join_values), family=('family', join_values), order=('order', join_values),
        candidate_read_count=('readid', 'nunique'), median_candidate_confidence=('species_bootstrap', 'median'),
        **{'class': ('class', join_values)})
    passed = candidates[candidates.species_bootstrap.ge(cutoff)]
    confident = passed.groupby(keys, dropna=False).agg(
        median_confidence=('species_bootstrap', 'median'), read_count=('readid', 'nunique'),
        raw_read_count=('raw_ont_read_id', 'nunique'),
        forward_primers_detected=('forward_PCR_primer', join_values),
        reverse_primers_detected=('reverse_PCR_primer', join_values))
    summary = all_species.join(confident).reset_index()
    for col in ['read_count', 'raw_read_count']:
        summary[col] = summary[col].fillna(0).astype(int)
    for col in ['forward_primers_detected', 'reverse_primers_detected']:
        summary[col] = summary[col].fillna('')
    summary['species_name'] = summary.species.map(lambda x: taxon_parts(x)[0])
    summary['species_taxid'] = summary.species.map(lambda x: taxon_parts(x)[1])
    for col in ['genus', 'family', 'order', 'class']:
        summary[col] = summary[col].map(lambda x: ';'.join(taxon_parts(v)[0] for v in x.split(';') if v))
    summary['contamination'] = summary.species.map(
        lambda label: 'contamination' if ('species', label) in flagged else 'nocontamination')
    return summary[columns].sort_values(['read_count', 'species'], ascending=[False, True])


def tag_status_rows(folder, tags, frame):
    """Tag outcomes counted over raw reads and over cleaned (classified) queries."""
    rows = []
    for direction in ('forward', 'reverse'):
        raw = tags[direction + '_status'].value_counts()
        cleaned = frame[direction + '_tag_status'].value_counts()
        for status in sorted(set(raw.index) | set(cleaned.index)):
            rows.append(dict(source_folder=folder, direction=direction, status=status,
                             raw_reads=int(raw.get(status, 0)), cleaned_queries=int(cleaned.get(status, 0))))
    return rows


def run_analysis(*, merged_dir, sintax_root, output_dir, folders=None, cutoff=0.8, prevalence=0.8,
                 suffix='sintaxSP', allow_missing=False):
    merged_dir, sintax_root, output_dir = dataset_path(merged_dir), Path(sintax_root), Path(output_dir)
    if folders is None:
        folders = [p.name for p in discover_samples([merged_dir])]
    if not 0 <= cutoff <= 1 or not 0 < prevalence <= 1:
        raise ValueError('Invalid confidence/prevalence threshold')
    if len(set(folders)) != len(folders):
        raise ValueError('Duplicate folders')
    if output_dir.exists():
        raise FileExistsError(f'Choose a new output directory: {output_dir}')
    inputs, manifest = {}, []
    for folder in folders:
        files = dict(sintax=sintax_root/folder/f'{folder}_{suffix}.tsv',
                     clean_reads=merged_dir/folder/'filtered/asinput.fastq.gz',
                     primer_tags=merged_dir/folder/'filtered/primer_tags.tsv.gz')
        row = dict(source_folder=folder, status='analysed' if files['sintax'].is_file() else 'missing_sintax')
        if row['status'] == 'analysed':
            for path in files.values():
                if not path.is_file():
                    raise FileNotFoundError(f'{path}; run Read processing for {folder}')
            inputs[folder] = files
        elif not allow_missing:
            raise FileNotFoundError(files['sintax'])
        for kind, path in files.items():
            row[kind + '_file'] = project_path(path)
        manifest.append(row)
    if not inputs:
        raise ValueError('No available SINTAX outputs')
    barcodes = {folder for folder in inputs if folder != 'unclassified'}
    output_dir.mkdir(parents=True)
    metadata = dict(created_utc=datetime.now(timezone.utc).isoformat(), cutoff=cutoff,
        prevalence_threshold=prevalence, denominator='distinct analysed barcode folders; unclassified excluded',
        analysed_barcode_count=len(barcodes), minimum_barcode_count=math.ceil(prevalence*len(barcodes)-1e-12),
        requested_folders=folders, analysed_folders=list(inputs),
        missing_folders=[r['source_folder'] for r in manifest if r['status'] != 'analysed'],
        merged_dir=project_path(merged_dir), sintax_root=project_path(sintax_root),
        input_files={folder: {kind: file_record(path) for kind, path in files.items()} for folder, files in inputs.items()},
        read_count_unit='unique cleaned query ID; adapter-split fragments are separate queries',
        confidence_source='reported rounded SINTAX bootstrap values, not calibrated probabilities')
    (output_dir/'analysis_parameters.json').write_text(json.dumps(metadata, indent=2)+'\n')
    write_tsv(pd.DataFrame(manifest), output_dir/'input_manifest.tsv')
    summaries, stats, tag_rows, presence = [], [], [], defaultdict(set)
    taxon_primers = {'forward': defaultdict(set), 'reverse': defaultdict(set)}
    with tempfile.TemporaryDirectory(prefix='sintax_analysis_') as temp:
        temp = Path(temp)
        for folder, files in inputs.items():
            print(f'Parsing {folder}: {files["sintax"]}', flush=True)
            tags = pd.read_csv(files['primer_tags'], sep='\t', dtype=str, keep_default_na=False,
                               usecols=['read_id', *TAG_FIELDS]).set_index('read_id')
            frame = sample_frame(files['sintax'], cutoff, tags)
            tag_rows += tag_status_rows(folder, tags, frame)
            for rank in RANK_NAMES:
                accepted = frame[frame['accepted_'+rank].ne('')]
                if folder != 'unclassified':
                    for label in accepted['accepted_'+rank].unique():
                        presence[(rank, label)].add(folder)
                for direction in ('forward', 'reverse'):
                    pairs = accepted[['accepted_'+rank, direction+'_PCR_primer']].drop_duplicates()
                    for label, primer in pairs.itertuples(index=False, name=None):
                        if primer:
                            taxon_primers[direction][(rank, label)].add(primer)
            frame.to_pickle(temp/f'{folder}.pkl')
            stats.append(dict(source_folder=folder, raw_read_count=len(tags),
                cleaned_query_count=len(frame), cleaned_raw_read_count=frame.raw_ont_read_id.nunique(),
                assigned_query_count=int(frame.assigned_label.ne('').sum()),
                species_query_count=int(frame.accepted_species.ne('').sum()),
                genus_or_species_lowest_query_count=int(frame.taxonomy_status.eq('taxapass').sum()),
                with_forward_primer_count=int(frame.forward_PCR_primer.ne('').sum()),
                with_reverse_primer_count=int(frame.reverse_PCR_primer.ne('').sum()),
                forward_primers_detected=join_values(frame.forward_PCR_primer),
                reverse_primers_detected=join_values(frame.reverse_PCR_primer)))
            del frame, tags
        prevalence_frame = prevalence_table(presence, barcodes, prevalence)
        for direction in ('forward', 'reverse'):
            prevalence_frame[direction+'_primers_detected'] = [
                ';'.join(sorted(taxon_primers[direction][(rank, label)]))
                for rank, label in prevalence_frame[['rank', 'taxon_label']].itertuples(index=False, name=None)]
        write_tsv(prevalence_frame, output_dir/'taxon_barcode_prevalence.tsv')
        flagged = set(prevalence_frame.loc[prevalence_frame.contamination.eq('contamination'),
                                           ['rank', 'taxon_label']].itertuples(index=False, name=None))
        print(f'Prevalence screen ({prevalence:.0%}): {len(barcodes)} barcodes, '
              f'{metadata["minimum_barcode_count"]} required; {len(flagged)} rank/taxon flags', flush=True)
        for folder in inputs:
            print(f'Writing read tables and species summary: {folder}', flush=True)
            frame = add_flags(pd.read_pickle(temp/f'{folder}.pkl'), flagged)
            folder_out = output_dir/folder
            write_tsv(frame, folder_out/f'{folder}_reads.tsv')
            summary = species_summary(frame, flagged, cutoff)
            summaries.append(summary)
            write_tsv(summary, folder_out/f'{folder}_species_summary.tsv')
            # Taxonomy per detected forward/reverse tag pair, including single-tag reads.
            primer_summary = (frame[frame.assigned_label.ne('') & (frame.forward_PCR_primer.ne('') | frame.reverse_PCR_primer.ne(''))]
                .groupby(['source_folder', 'forward_PCR_primer', 'reverse_PCR_primer',
                          'assigned_rank', 'assigned_label', 'contamination'], dropna=False)
                .agg(read_count=('readid', 'nunique'), raw_read_count=('raw_ont_read_id', 'nunique'),
                     median_confidence=('assigned_confidence', 'median')).reset_index())
            write_tsv(primer_summary, folder_out/f'{folder}_primer_taxonomy.tsv')
            (temp/f'{folder}.pkl').unlink()
            del frame
    write_tsv(pd.concat(summaries, ignore_index=True), output_dir/'all_species_summary.tsv')
    write_tsv(pd.DataFrame(stats), output_dir/'barcode_analysis_summary.tsv')
    write_tsv(pd.DataFrame(tag_rows), output_dir/'primer_tag_summary.tsv')
    write_report(output_dir)
    describe_headers(output_dir)
    (output_dir/'analysis.complete').write_text(datetime.now(timezone.utc).isoformat()+'\n')
    print(f'Analysis complete: {output_dir}', flush=True)
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('dataset', type=Path, help='Dataset parent, e.g. sequencing2 or data/sequencing2')
    parser.add_argument('--output-dir', type=Path, required=True, help='New directory for the analysis tables')
    parser.add_argument('--sintax-root', type=Path, help='Default: data/sintax/DATASET')
    parser.add_argument('--folders', nargs='+', help='Optional sample-name subset; default: all samples in the dataset')
    parser.add_argument('--cutoff', type=float, default=0.8)
    parser.add_argument('--prevalence', type=float, default=0.8)
    parser.add_argument('--suffix', default='sintaxSP')
    parser.add_argument('--allow-missing', action='store_true', help='Record samples without SINTAX output instead of failing')
    args = parser.parse_args()
    dataset = dataset_path(args.dataset)
    sintax_root = args.sintax_root or Path(__file__).resolve().parents[1]/'data/sintax'/dataset.name
    run_analysis(merged_dir=dataset, sintax_root=sintax_root, output_dir=args.output_dir, folders=args.folders,
                 cutoff=args.cutoff, prevalence=args.prevalence, suffix=args.suffix, allow_missing=args.allow_missing)


if __name__ == '__main__':
    main()
