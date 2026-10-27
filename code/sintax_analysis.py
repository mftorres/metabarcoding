#!/usr/bin/env python3
"""Read-level SINTAX analysis and PCR-plate prevalence, without classification."""
import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import tempfile

import pandas as pd
from sintax_parser import RANKS, COLUMNS, parse_output, taxon_parts

RANK_NAMES = list(RANKS.values())
PRIMER_COLUMNS = ['merged_folder', 'barcode_folder', 'seqID', 'forward_patternName', 'reverse_patternName']


def join_values(values):
    return ';'.join(sorted({name for v in values if pd.notna(v) for name in str(v).split(';') if name}))


def write_tsv(frame, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep='\t', index=False, header=True, mode='x', na_rep='')


def primer_name_map(mft_fasta, pcr_fasta):
    def records(path):
        result, name, sequence = {}, None, []
        with Path(path).open() as handle:
            for line in handle:
                if line.startswith('>'):
                    if name: result[name] = ''.join(sequence).upper()
                    name, sequence = line[1:].split()[0], []
                else: sequence.append(line.strip())
        if name: result[name] = ''.join(sequence).upper()
        return result
    mft, pcr = records(mft_fasta), records(pcr_fasta)
    by_sequence = defaultdict(list)
    for name, sequence in pcr.items(): by_sequence[sequence].append(name)
    return {name: ';'.join(sorted(by_sequence[seq])) for name, seq in mft.items() if seq in by_sequence}


def split_primer_table(primer_path, dataset, barcodes, temp_root):
    """Read a large best-hit table once, retaining only fields used in this analysis."""
    paths = {}
    counts = dict.fromkeys(barcodes, 0)
    for chunk in pd.read_csv(primer_path, sep='\t', dtype=str, usecols=PRIMER_COLUMNS,
                             keep_default_na=False, chunksize=100000):
        chunk = chunk[chunk.merged_folder.eq(dataset) & chunk.barcode_folder.isin(barcodes)]
        for barcode, group in chunk.groupby('barcode_folder', sort=False):
            path = temp_root / f'{barcode}_primers.tsv'
            group.to_csv(path, sep='\t', index=False, mode='a', header=barcode not in paths)
            paths[barcode] = path
            counts[barcode] += len(group)
    for barcode in barcodes:
        if barcode not in paths:
            raise ValueError(f'No primer tracking rows for {dataset}/{barcode}; cannot label missing tracking as absent primers')
    return paths, counts


def sample_frame(path, cutoff, primer_path, mapping, plate, dataset):
    frame = pd.DataFrame.from_records(parse_output(path, cutoff=cutoff), columns=COLUMNS)
    if frame.readid.duplicated().any():
        raise ValueError(f'Duplicate query IDs in {path}')
    for col in ['assigned_confidence', 'taxonomy_lowestboot'] + [r+'_bootstrap' for r in RANK_NAMES]:
        frame[col] = pd.to_numeric(frame[col], errors='coerce')
    frame['PCR_plate'] = plate
    frame['dataset'] = dataset
    primers = pd.read_csv(primer_path, sep='\t', dtype=str, keep_default_na=False)
    if primers.seqID.duplicated().any(): raise ValueError(f'Duplicate primer IDs in {primer_path}')
    primers = primers.set_index('seqID')
    # These primer results were generated on raw reads. Match by their parent ID.
    frame['primer_tracking_record_found'] = frame.raw_ont_read_id.isin(primers.index)
    for direction in ('forward', 'reverse'):
        column = direction+'_patternName'
        original_names = frame.raw_ont_read_id.map(primers[column]).fillna('')
        canonical_names = {name for value in mapping.values() for name in value.split(';')}
        unknown = set(original_names[original_names.ne('')]) - set(mapping) - canonical_names
        if unknown:
            raise ValueError(f'Unmapped {direction} primer names: {sorted(unknown)}')
        frame[direction+'_PCR_primer'] = original_names.map(mapping).fillna('')
        already_mapped = original_names.isin(canonical_names)
        frame.loc[already_mapped, direction+'_PCR_primer'] = original_names[already_mapped]
    return frame


def prevalence_table(presence, plates, fraction):
    minimum = math.ceil(fraction * len(plates) - 1e-12)
    rows = []
    for (rank, label), observed in sorted(presence.items()):
        name, taxid = taxon_parts(label)
        rows.append(dict(rank=rank, taxon_label=label, taxon_name=name, taxid=taxid,
                         plate_count=len(observed), analysed_plate_count=len(plates),
                         plate_fraction=len(observed)/len(plates), minimum_plate_count=minimum,
                         prevalence_threshold=fraction, PCR_plates=';'.join(sorted(observed)),
                         contamination='contamination' if len(observed) >= minimum else 'nocontamination'))
    return pd.DataFrame(rows, columns=['rank','taxon_label','taxon_name','taxid','plate_count',
        'analysed_plate_count','plate_fraction','minimum_plate_count','prevalence_threshold',
        'PCR_plates','contamination'])


def add_flags(frame, flagged):
    frame['contamination'] = [
        'contamination' if (rank, label) in flagged else 'nocontamination'
        for rank, label in zip(frame.assigned_rank, frame.assigned_label)]
    for rank in RANK_NAMES:
        frame[rank+'_contamination'] = frame['accepted_'+rank].map(
            lambda label: 'contamination' if (rank, label) in flagged else 'nocontamination')
    return frame


def species_summary(frame, flagged, cutoff):
    """One row per candidate species; statistics distinguish confident and candidate calls."""
    candidates = frame[frame.species.ne('')].copy()
    keys = ['source_folder', 'source_name', 'source_file', 'PCR_plate', 'species']
    columns = keys + ['species_name','species_taxid','genus','family','order','class',
                     'contamination','median_confidence','read_count','raw_read_count',
                     'forward_primers_detected','reverse_primers_detected',
                     'candidate_read_count','median_candidate_confidence']
    if candidates.empty: return pd.DataFrame(columns=columns)
    all_species = candidates.groupby(keys, dropna=False).agg(
        genus=('genus',join_values), family=('family',join_values), order=('order',join_values),
        candidate_read_count=('readid','nunique'), median_candidate_confidence=('species_bootstrap','median'),
        **{'class': ('class', join_values)})
    passed = candidates[candidates.species_bootstrap.ge(cutoff)]
    confident = passed.groupby(keys, dropna=False).agg(
        median_confidence=('species_bootstrap','median'), read_count=('readid','nunique'),
        raw_read_count=('raw_ont_read_id','nunique'),
        forward_primers_detected=('forward_PCR_primer',join_values),
        reverse_primers_detected=('reverse_PCR_primer',join_values))
    summary = all_species.join(confident).reset_index()
    for col in ['read_count','raw_read_count']:
        summary[col] = summary[col].fillna(0).astype(int)
    for col in ['forward_primers_detected','reverse_primers_detected']:
        summary[col] = summary[col].fillna('')
    summary['species_name'] = summary.species.map(lambda x: taxon_parts(x)[0])
    summary['species_taxid'] = summary.species.map(lambda x: taxon_parts(x)[1])
    for col in ['genus','family','order','class']:
        summary[col] = summary[col].map(lambda x: ';'.join(taxon_parts(v)[0] for v in x.split(';') if v))
    summary['contamination'] = summary.species.map(
        lambda label: 'contamination' if ('species',label) in flagged else 'nocontamination')
    return summary[columns].sort_values(['read_count','species'],ascending=[False,True])


def run_analysis(*, project_root, merged_dir, sintax_root, output_dir, primer_hits,
                 plate_map, folders, cutoff=0.8, prevalence=0.8, suffix='sintaxSP', allow_missing=False):
    project_root, merged_dir, sintax_root, output_dir = map(Path, (project_root, merged_dir, sintax_root, output_dir))
    if not 0 <= cutoff <= 1 or not 0 < prevalence <= 1: raise ValueError('Invalid confidence/prevalence threshold')
    if len(set(folders)) != len(folders): raise ValueError('Duplicate folders')
    if output_dir.exists(): raise FileExistsError(f'Choose a new output directory: {output_dir}')
    plate_frame = pd.read_csv(plate_map, sep='\t', dtype=str, keep_default_na=False)
    plate_frame = plate_frame[['sequencing_barcode','PCR_plate']].drop_duplicates()
    if plate_frame.sequencing_barcode.duplicated().any(): raise ValueError('Ambiguous PCR_plate mapping')
    plate_lookup = dict(zip(plate_frame.sequencing_barcode,plate_frame.PCR_plate))
    inputs, manifest = {}, []
    for folder in folders:
        source = sintax_root / folder / f'{folder}_{suffix}.tsv'
        clean = merged_dir / folder / 'filtered/asinput.fastq.gz'
        status = 'analysed' if source.is_file() else 'missing_sintax'
        if source.is_file():
            if not clean.is_file(): raise FileNotFoundError(clean)
            if not plate_lookup.get(folder): raise ValueError(f'No PCR plate for {folder}')
            inputs[folder] = source
        elif not allow_missing: raise FileNotFoundError(source)
        manifest.append(dict(source_folder=folder, PCR_plate=plate_lookup.get(folder,''),
            status=status, sintax_file=str(source.resolve()), clean_reads_file=str(clean.resolve()),
            sintax_bytes=source.stat().st_size if source.exists() else '',
            sintax_mtime_ns=source.stat().st_mtime_ns if source.exists() else ''))
    if not inputs: raise ValueError('No available SINTAX outputs')
    plates = {plate_lookup[folder] for folder in inputs}
    mapping = primer_name_map(project_root/'data/saproxdna_primertags_seq5_2026_1.fasta',
                              project_root/'data/PCRsaproxdna_primertags_seq5_2026_1.fasta')
    output_dir.mkdir(parents=True)
    metadata = dict(created_utc=datetime.now(timezone.utc).isoformat(), cutoff=cutoff,
        prevalence_threshold=prevalence, denominator='distinct PCR_plate codes of analysed folders',
        analysed_plate_count=len(plates), minimum_plate_count=math.ceil(prevalence*len(plates)-1e-12),
        requested_folders=folders, analysed_folders=list(inputs),
        missing_folders=[r['source_folder'] for r in manifest if r['status']!='analysed'],
        primer_hits=str(Path(primer_hits).resolve()), plate_map=str(Path(plate_map).resolve()),
        merged_dir=str(merged_dir.resolve()), sintax_root=str(sintax_root.resolve()),
        read_count_unit='unique cleaned query ID; adapter fragments are separate queries',
        median_confidence_scope='species confidence among queries with species confidence >= cutoff',
        confidence_source='reported rounded SINTAX bootstrap values, not calibrated probabilities')
    (output_dir/'analysis_parameters.json').write_text(json.dumps(metadata,indent=2)+'\n')
    write_tsv(pd.DataFrame(manifest),output_dir/'input_manifest.tsv')
    (output_dir/'primer_name_mapping.json').write_text(json.dumps(mapping, indent=2)+'\n')
    summaries, stats, presence = [], [], defaultdict(set)
    taxon_primers = {'forward': defaultdict(set), 'reverse': defaultdict(set)}
    with tempfile.TemporaryDirectory(prefix='sintax_analysis_') as temp:
        temp = Path(temp)
        print('Partitioning primer best hits (one pass)...',flush=True)
        primer_paths, primer_counts = split_primer_table(primer_hits,merged_dir.name,list(inputs),temp)
        for folder, path in inputs.items():
            print(f'Parsing {folder}: {path}',flush=True)
            frame = sample_frame(path,cutoff,primer_paths[folder],mapping,plate_lookup[folder],merged_dir.name)
            for rank in RANK_NAMES:
                accepted = frame[frame['accepted_'+rank].ne('')]
                for label in accepted['accepted_'+rank].unique():
                    presence[(rank,label)].add(plate_lookup[folder])
                for direction in ('forward', 'reverse'):
                    for label, primer in accepted[['accepted_'+rank, direction+'_PCR_primer']].drop_duplicates().itertuples(index=False, name=None):
                        taxon_primers[direction][(rank,label)].update(name for name in primer.split(';') if name)
            frame.to_pickle(temp/f'{folder}.pkl')
            stats.append(dict(source_folder=folder,PCR_plate=plate_lookup[folder],
                cleaned_query_count=len(frame),original_read_count=frame.raw_ont_read_id.nunique(),
                assigned_query_count=int(frame.assigned_label.ne('').sum()),
                species_query_count=int(frame.accepted_species.ne('').sum()),
                genus_or_species_lowest_query_count=int(frame.taxonomy_status.eq('taxapass').sum()),
                with_forward_primer_count=int(frame.forward_PCR_primer.ne('').sum()),
                with_reverse_primer_count=int(frame.reverse_PCR_primer.ne('').sum()),
                forward_primers_detected=join_values(frame.forward_PCR_primer),
                reverse_primers_detected=join_values(frame.reverse_PCR_primer),
                primer_tracking_rows=primer_counts[folder],
                parent_reads_with_primer_tracking=int(frame.primer_tracking_record_found.sum())))
            del frame
        prevalence_frame = prevalence_table(presence,plates,prevalence)
        for direction in ('forward', 'reverse'):
            prevalence_frame[direction+'_primers_detected'] = [
                ';'.join(sorted(taxon_primers[direction][(rank,label)]))
                for rank,label in prevalence_frame[['rank','taxon_label']].itertuples(index=False,name=None)]
        write_tsv(prevalence_frame,output_dir/'taxon_plate_prevalence.tsv')
        flagged = set(prevalence_frame.loc[prevalence_frame.contamination.eq('contamination'), ['rank','taxon_label']].itertuples(index=False,name=None))
        print(f'80% screen: {len(plates)} plates, {metadata["minimum_plate_count"]} required; {len(flagged)} rank/taxon flags',flush=True)
        combined = output_dir/'all_reads.tsv'
        with combined.open('x',encoding='utf-8',newline='') as combined_handle:
            for i,folder in enumerate(inputs):
                print(f'Writing read tables and species summary: {folder}',flush=True)
                frame = add_flags(pd.read_pickle(temp/f'{folder}.pkl'),flagged)
                folder_out = output_dir/folder
                write_tsv(frame,folder_out/f'{folder}_reads.tsv')
                frame.to_csv(combined_handle,sep='\t',index=False,header=i==0,na_rep='')
                summary = species_summary(frame,flagged,cutoff)
                summaries.append(summary)
                write_tsv(summary,folder_out/f'{folder}_species_summary.tsv')
                # Per-primer taxonomy evidence retains folder and plate provenance.
                primer_summary = (frame[frame.assigned_label.ne('') & (frame.forward_PCR_primer.ne('') | frame.reverse_PCR_primer.ne(''))]
                    .groupby(['source_folder','PCR_plate','forward_PCR_primer','reverse_PCR_primer',
                              'assigned_rank','assigned_label','contamination'],dropna=False)
                    .agg(read_count=('readid','nunique'),raw_read_count=('raw_ont_read_id','nunique'),
                         median_confidence=('assigned_confidence','median')).reset_index())
                write_tsv(primer_summary,folder_out/f'{folder}_primer_taxonomy.tsv')
                (temp/f'{folder}.pkl').unlink()
                del frame
    write_tsv(pd.concat(summaries,ignore_index=True),output_dir/'all_species_summary.tsv')
    write_tsv(pd.DataFrame(stats),output_dir/'barcode_analysis_summary.tsv')
    (output_dir/'analysis.complete').write_text(datetime.now(timezone.utc).isoformat()+'\n')
    print(f'Analysis complete: {output_dir}',flush=True)
    return metadata


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--merged-dir',type=Path,default=Path('data/mergedseq5'))
    parser.add_argument('--sintax-root',type=Path,default=Path('data/sintax'))
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--primer-hits',type=Path,required=True)
    parser.add_argument('--plate-map',type=Path,default=Path('saproxDNA/species_info_seq5/seq5_pcr_to_pools.tsv'))
    parser.add_argument('--folders',nargs='+',default=[f'barcode{i:02d}' for i in range(1,25)])
    parser.add_argument('--cutoff',type=float,default=0.8)
    parser.add_argument('--prevalence',type=float,default=0.8)
    parser.add_argument('--suffix',default='sintaxSP')
    parser.add_argument('--allow-missing',action='store_true')
    args=vars(parser.parse_args())
    run_analysis(**args)

if __name__=='__main__': main()
