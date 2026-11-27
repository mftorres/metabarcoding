#!/usr/bin/env python3
"""Module 1 analysis tables: taxonomy of every cluster, species and specimen summaries,
and the contamination screen, for one dataset.

Inputs per sample folder: the clusters and their members (cluster_specimens.py), the
tag table (primer_tags.py) and the SINTAX classification of the clusters. Every
count of reads is the number of cleaned reads in the clusters concerned.
"""
import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import tempfile

import pandas as pd

from analysis_report import write_report
from sample_folders import dataset_path, discover_samples, pipeline_setting, project_path, PROJECT_ROOT
from seqio import fasta, raw_read_id
from sintax_parser import RANKS, parse_output, taxon_parts

RANK_NAMES = list(RANKS.values())
CLUSTER_COLUMNS = (['source_folder', 'cluster_id', 'cluster_type', 'read_count', 'forward_tag', 'reverse_tags',
                    'assigned_rank', 'assigned_label', 'assigned_name', 'assigned_taxid', 'assigned_confidence',
                    'contamination'] + [column for rank in RANK_NAMES for column in (rank, rank + '_confidence')]
                   + ['strand', 'sintax_taxonomy'])


def write_tsv(frame, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep='\t', index=False, mode='x', na_rep='')


def file_record(path):
    return dict(path=project_path(path), bytes=Path(path).stat().st_size)


def join_values(values):
    """Unique non-empty values as one semicolon-separated string."""
    return ';'.join(sorted({value for value in values if value}))


def cluster_table(folder, files, cutoff):
    """One row per cluster: read count, tags and the ranks SINTAX called with confidence >= cutoff."""
    members = pd.read_csv(files['cluster_members'], sep='\t', dtype=str, keep_default_na=False)
    tags = pd.read_csv(files['primer_tags'], sep='\t', dtype=str, keep_default_na=False,
                       usecols=['read_id', 'forward_status', 'reverse_primer', 'reverse_status']).set_index('read_id')
    members['raw_read_id'] = members.read_id.map(raw_read_id)
    members['reverse_tag'] = members.raw_read_id.map(tags.reverse_primer).fillna('')
    clusters = members.groupby('cluster_id').agg(read_count=('read_id', 'size'), forward_tag=('forward_tag', 'first'),
                                                 reverse_tags=('reverse_tag', join_values))
    sequences = {name for name, _ in fasta(files['clusters'])}
    if sequences != set(clusters.index):
        raise ValueError(f'{folder}: clusters.fasta and cluster_members.tsv.gz do not match; rerun the clustering')

    # SINTAX skips sequences outside its length limits: those clusters stay unassigned.
    rows = {}
    for name, strand, taxonomy, ranks in parse_output(files['sintax']):
        if name not in sequences:
            raise ValueError(f'{files["sintax"]}: {name} is not one of the clusters; '
                             'the classification is older than the clusters, rerun it (FORCE_SINTAX=1)')
        row = dict(strand=strand, sintax_taxonomy=taxonomy)
        for rank in RANK_NAMES:
            label, confidence = ranks.get(rank, ('', None))
            row[rank] = label if confidence is not None and confidence >= cutoff else ''
            row[rank + '_confidence'] = confidence
            # The lowest accepted rank is the cluster's identification.
            if row[rank]:
                name_only, taxid = taxon_parts(label)
                row.update(assigned_rank=rank, assigned_label=label, assigned_name=name_only,
                           assigned_taxid=taxid, assigned_confidence=confidence)
        rows[name] = row

    frame = pd.DataFrame.from_dict(rows, orient='index').reindex(clusters.index)
    frame = clusters.join(frame).rename_axis('cluster_id').reset_index()
    frame['source_folder'] = folder
    frame['cluster_type'] = ['consensus' if count >= 2 else 'single_read' for count in frame.read_count]
    frame['contamination'] = ''
    frame = frame.reindex(columns=CLUSTER_COLUMNS)
    text_columns = [column for column in CLUSTER_COLUMNS if column not in ('read_count',) and 'confidence' not in column]
    frame[text_columns] = frame[text_columns].fillna('')
    return frame, members, tags


def tag_status_rows(folder, tags, members):
    """Tag outcomes counted over raw reads and over cleaned reads."""
    rows = []
    cleaned = members.raw_read_id.map(tags.forward_status), members.raw_read_id.map(tags.reverse_status)
    for direction, cleaned_status in zip(('forward', 'reverse'), cleaned):
        raw_counts = tags[direction + '_status'].value_counts()
        cleaned_counts = cleaned_status.value_counts()
        for status in sorted(set(raw_counts.index) | set(cleaned_counts.index)):
            rows.append(dict(source_folder=folder, direction=direction, status=status,
                             raw_reads=int(raw_counts.get(status, 0)), cleaned_reads=int(cleaned_counts.get(status, 0))))
    return rows


def prevalence_table(presence, barcodes, fraction):
    """Barcode folders in which each accepted taxon occurs; frequent taxa are flagged."""
    minimum = math.ceil(fraction * len(barcodes) - 1e-12)
    rows = []
    for (rank, label), found in sorted(presence.items()):
        name, taxid = taxon_parts(label)
        rows.append(dict(rank=rank, taxon_label=label, taxon_name=name, taxid=taxid,
                         barcode_count=len(found), analysed_barcode_count=len(barcodes),
                         barcode_fraction=round(len(found) / len(barcodes), 4), minimum_barcode_count=minimum,
                         barcodes=';'.join(sorted(found)),
                         contamination='contamination' if len(found) >= minimum else 'nocontamination'))
    return pd.DataFrame(rows, columns=['rank', 'taxon_label', 'taxon_name', 'taxid', 'barcode_count',
                                       'analysed_barcode_count', 'barcode_fraction', 'minimum_barcode_count',
                                       'barcodes', 'contamination'])


def species_summary(frame):
    """One row per species with confident clusters in a folder."""
    species = frame[frame.species.ne('')]
    summary = species.groupby(['source_folder', 'species']).agg(
        read_count=('read_count', 'sum'), cluster_count=('cluster_id', 'size'),
        specimen_count=('forward_tag', lambda tags: len(set(tags) - {'none'})),
        forward_tags=('forward_tag', join_values), median_confidence=('species_confidence', 'median'),
        genus=('genus', join_values), family=('family', join_values), order=('order', join_values),
        contamination=('contamination', join_values), **{'class': ('class', join_values)}).reset_index()
    summary.insert(2, 'species_name', [taxon_parts(label)[0] for label in summary.species])
    summary.insert(3, 'species_taxid', [taxon_parts(label)[1] for label in summary.species])
    for rank in ('genus', 'family', 'order', 'class'):
        summary[rank] = [';'.join(taxon_parts(label)[0] for label in value.split(';') if label) for value in summary[rank]]
    return summary.sort_values(['read_count', 'species'], ascending=[False, True])


def specimen_taxa(frame):
    """Identifications within each specimen (barcode folder + forward tag), largest first."""
    taxa = frame.assign(assigned_label=frame.assigned_label.replace('', 'unassigned'))
    table = taxa.groupby(['source_folder', 'forward_tag', 'assigned_rank', 'assigned_label', 'contamination']).agg(
        read_count=('read_count', 'sum'), cluster_count=('cluster_id', 'size')).reset_index()
    totals = table.groupby(['source_folder', 'forward_tag']).read_count.transform('sum')
    table['fraction_of_specimen_reads'] = (table.read_count / totals).round(4)
    largest = table.groupby(['source_folder', 'forward_tag']).read_count.transform('max')
    table['largest_taxon'] = table.read_count.eq(largest)
    return table.sort_values(['source_folder', 'forward_tag', 'read_count'], ascending=[True, True, False])


def run_analysis(*, dataset, sintax_root, output_dir, folders=None, cutoff=None, prevalence=None, suffix=None,
                 allow_missing=False):
    dataset, sintax_root, output_dir = dataset_path(dataset), Path(sintax_root), Path(output_dir)
    cutoff = float(cutoff if cutoff is not None else pipeline_setting('ANALYSIS_CUTOFF'))
    prevalence = float(prevalence if prevalence is not None else pipeline_setting('PREVALENCE'))
    suffix = suffix or pipeline_setting('SINTAX_SUFFIX')
    if folders is None:
        folders = [path.name for path in discover_samples([dataset])]
    if not 0 <= cutoff <= 1 or not 0 < prevalence <= 1:
        raise ValueError('The confidence cutoff and prevalence must be between 0 and 1')
    if output_dir.exists():
        raise FileExistsError(f'Choose a new output directory: {output_dir}')

    # Find the input files of every folder.
    inputs, missing = {}, []
    for folder in folders:
        filtered = dataset/folder/'filtered'
        files = dict(sintax=sintax_root/folder/f'{folder}_{suffix}.tsv', clusters=filtered/'clusters.fasta',
                     cluster_members=filtered/'cluster_members.tsv.gz', primer_tags=filtered/'primer_tags.tsv.gz')
        if not files['sintax'].is_file():
            if not allow_missing:
                raise FileNotFoundError(f'{files["sintax"]}; run Module 1 for {folder}')
            missing.append(folder)
            continue
        for path in files.values():
            if not path.is_file():
                raise FileNotFoundError(f'{path}; run Module 1 for {folder}')
        inputs[folder] = files
    if not inputs:
        raise ValueError('No folder has a SINTAX classification')
    barcodes = [folder for folder in inputs if folder != 'unclassified']
    output_dir.mkdir(parents=True)
    metadata = dict(created_utc=datetime.now(timezone.utc).isoformat(), dataset=project_path(dataset),
                    sintax_root=project_path(sintax_root), sintax_suffix=suffix,
                    cluster_identity=pipeline_setting('CLUSTER_IDENTITY'), confidence_cutoff=cutoff,
                    prevalence_threshold=prevalence, analysed_barcode_count=len(barcodes),
                    minimum_barcode_count=math.ceil(prevalence * len(barcodes) - 1e-12),
                    analysed_folders=list(inputs), missing_folders=missing,
                    input_files={folder: {kind: file_record(path) for kind, path in files.items()}
                                 for folder, files in inputs.items()})
    (output_dir/'analysis_parameters.json').write_text(json.dumps(metadata, indent=2) + '\n')

    # Pass 1: cluster tables, read counts and the barcodes in which each accepted taxon occurs.
    # Tables wait in a temporary folder so that only one sample is held in memory.
    stats, tag_rows, presence = [], [], {}
    with tempfile.TemporaryDirectory(prefix='sintax_analysis_') as temp:
        temp = Path(temp)
        for folder, files in inputs.items():
            print(f'Reading {folder}', flush=True)
            frame, members, tags = cluster_table(folder, files, cutoff)
            tag_rows += tag_status_rows(folder, tags, members)
            if folder != 'unclassified':
                for rank in RANK_NAMES:
                    for label in frame.loc[frame[rank].ne(''), rank].unique():
                        presence.setdefault((rank, label), set()).add(folder)
            reads = frame.read_count
            stats.append(dict(source_folder=folder, raw_reads=len(tags), cleaned_reads=int(reads.sum()),
                              reads_without_forward_tag=int(reads[frame.forward_tag.eq('none')].sum()),
                              specimens=frame.loc[frame.forward_tag.ne('none'), 'forward_tag'].nunique(),
                              clusters=len(frame), consensus_sequences=int(frame.cluster_type.eq('consensus').sum()),
                              single_reads=int(frame.cluster_type.eq('single_read').sum()),
                              reads_assigned=int(reads[frame.assigned_label.ne('')].sum()),
                              reads_to_genus=int(reads[frame.genus.ne('')].sum()),
                              reads_to_species=int(reads[frame.species.ne('')].sum()),
                              species=frame.loc[frame.species.ne(''), 'species'].nunique()))
            frame.to_pickle(temp/f'{folder}.pkl')

        prevalence_frame = prevalence_table(presence, barcodes, prevalence)
        write_tsv(prevalence_frame, output_dir/'taxon_barcode_prevalence.tsv')
        flagged = set(prevalence_frame.loc[prevalence_frame.contamination.eq('contamination'), ['rank', 'taxon_label']]
                      .itertuples(index=False, name=None))
        print(f'Contamination screen: taxa in at least {metadata["minimum_barcode_count"]} of {len(barcodes)} '
              f'barcode folders are flagged ({len(flagged)} rank/taxon pairs)', flush=True)

        # Pass 2: flag clusters whose identification is a frequent taxon, then write the tables.
        summaries, specimens = [], []
        for folder, row in zip(inputs, stats):
            frame = pd.read_pickle(temp/f'{folder}.pkl')
            frame['contamination'] = ['contamination' if (rank, label) in flagged else 'nocontamination'
                                      for rank, label in zip(frame.assigned_rank, frame.assigned_label)]
            flagged_rows = frame.contamination.eq('contamination')
            row['reads_flagged_contamination'] = int(frame.read_count[flagged_rows].sum())
            write_tsv(frame, output_dir/folder/f'{folder}_clusters.tsv')
            summaries.append(species_summary(frame))
            specimens.append(specimen_taxa(frame))
    write_tsv(pd.concat(summaries, ignore_index=True), output_dir/'all_species_summary.tsv')
    write_tsv(pd.concat(specimens, ignore_index=True), output_dir/'specimen_taxa.tsv')
    write_tsv(pd.DataFrame(stats), output_dir/'barcode_analysis_summary.tsv')
    write_tsv(pd.DataFrame(tag_rows), output_dir/'primer_tag_summary.tsv')
    write_report(output_dir)
    (output_dir/'analysis.complete').write_text(datetime.now(timezone.utc).isoformat() + '\n')
    print(f'Analysis complete: {output_dir}', flush=True)
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('dataset', help='Dataset folder, e.g. mergedseq7 or data/mergedseq7')
    parser.add_argument('--output-dir', type=Path, required=True, help='New folder for the analysis tables')
    parser.add_argument('--sintax-root', type=Path, help='Classification folder; default data/sintax/DATASET')
    parser.add_argument('--folders', nargs='+', help='Only these sample folders; default all')
    parser.add_argument('--cutoff', type=float, help='Default ANALYSIS_CUTOFF in config/pipeline.env')
    parser.add_argument('--prevalence', type=float, help='Default PREVALENCE in config/pipeline.env')
    parser.add_argument('--suffix', help='Classification file suffix; default SINTAX_SUFFIX in config/pipeline.env')
    parser.add_argument('--allow-missing', action='store_true', help='Skip folders without a classification')
    args = parser.parse_args()
    dataset = dataset_path(args.dataset)
    run_analysis(dataset=dataset, sintax_root=args.sintax_root or PROJECT_ROOT/'data/sintax'/dataset.name,
                 output_dir=args.output_dir, folders=args.folders, cutoff=args.cutoff, prevalence=args.prevalence,
                 suffix=args.suffix, allow_missing=args.allow_missing)


if __name__ == '__main__':
    main()
