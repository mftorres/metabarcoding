#!/usr/bin/env python3
"""Select SINTAX taxa with queried-rank confidence >=0.80; annotate with the lowest trusted taxon."""
import argparse
import csv
import difflib
import gzip
import json
import os
from pathlib import Path
import re
import tempfile
from sequence_headers import header_component, export_header, primer_annotations
from sintax_parser import RANKS, parse_output, taxon_parts

TAXON_SELECTION_POLICY = "queried-rank-confidence-ge-0.80-v2"


def normalize(value):
    return re.sub(r'[\s_]+',' ',value.strip()).casefold()


def fastq_records(path):
    opener = gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt') as handle:
        while True:
            header=handle.readline()
            if not header: break
            sequence=handle.readline().rstrip('\r\n')
            plus=handle.readline()
            quality=handle.readline().rstrip('\r\n')
            if not header.startswith('@') or not plus.startswith('+') or len(sequence)!=len(quality):
                raise ValueError(f'Invalid FASTQ record in {path}: {header[:100]!r}')
            yield header[1:].split()[0],sequence


def export_reads(*,taxa,folders,analysis_dir,output,merged_dir=None,sintax_root=None,suffix='sintaxSP'):
    analysis_dir,output=Path(analysis_dir),Path(output)
    if output.exists(): raise FileExistsError(output)
    if not (analysis_dir/'analysis.complete').is_file(): raise ValueError('Analysis is incomplete')
    params=json.loads((analysis_dir/'analysis_parameters.json').read_text())
    merged_dir=Path(merged_dir or params['merged_dir'])
    sintax_root=Path(sintax_root or params['sintax_root'])
    cutoff=params['cutoff']
    if cutoff != 0.8:
        raise ValueError('Headers require a completed analysis with cutoff 0.80')
    with (analysis_dir/'taxon_plate_prevalence.tsv').open() as handle:
        flagged={(r['rank'],r['taxon_label']) for r in csv.DictReader(handle,delimiter='\t')
                 if r['contamination']=='contamination'}
    with (analysis_dir/'input_manifest.tsv').open() as handle:
        manifest={r['source_folder']:r for r in csv.DictReader(handle,delimiter='\t')}
    names=[Path(f).name for f in folders]
    if len(set(names))!=len(names): raise ValueError('Duplicate folder names')
    targets={normalize(t) for t in taxa}
    seen_taxa=set()
    total=0
    output.parent.mkdir(parents=True,exist_ok=True)
    # Build privately and publish only after all matched cleaned reads are found.
    with tempfile.TemporaryDirectory(prefix='.sintax_fasta_',dir=output.parent) as temp:
        temporary=Path(temp)/'reads.fasta'
        with temporary.open('w') as writer:
            for folder in names:
                if folder not in params['analysed_folders']:
                    raise ValueError(f'{folder} was not analysed; contamination context is unavailable')
                source=sintax_root/folder/f'{folder}_{suffix}.tsv'
                clean=merged_dir/folder/'filtered/asinput.fastq.gz'
                if str(source.resolve())!=manifest[folder]['sintax_file'] or str(source.stat().st_size)!=manifest[folder]['sintax_bytes'] or str(source.stat().st_mtime_ns)!=manifest[folder]['sintax_mtime_ns']:
                    raise ValueError(f'SINTAX file changed since analysis: {source}; rerun analyses first')
                selected={}
                for row in parse_output(source,cutoff=cutoff):
                    matches=False
                    for rank in RANKS.values():
                        label=row[rank]
                        if not label: continue
                        name,taxid=taxon_parts(label)
                        seen_taxa.add(normalize(name))
                        confidence = row[rank + '_bootstrap']
                        if confidence != '' and confidence >= cutoff and targets.intersection({normalize(name),normalize(label),taxid}):
                            matches=True
                    if matches:
                        if row['readid'] in selected: raise ValueError(f'Duplicate query ID: {row["readid"]}')
                        flag='contamination' if (row['assigned_rank'],row['assigned_label']) in flagged else 'nocontamination'
                        selected[row['readid']]=(row, flag)
                primers = primer_annotations(analysis_dir, folder, {q: row['raw_ont_read_id'] for q, (row, _) in selected.items()})
                selected = {q: export_header(folder, *primers[q], row['assigned_label'], flag, q)
                            for q, (row, flag) in selected.items()}
                wanted=len(selected)
                for readid,sequence in fastq_records(clean):
                    label=selected.pop(readid,None)
                    if label is not None:
                        writer.write(f'>{label}\n')
                        for start in range(0,len(sequence),80): writer.write(sequence[start:start+80]+'\n')
                        total+=1
                if selected:
                    raise ValueError(f'{folder}: {len(selected)} selected query IDs absent from cleaned reads; examples: {list(selected)[:3]}')
                print(f'{folder}: exported {wanted:,} cleaned reads',flush=True)
        if not total:
            suggestions={t:difflib.get_close_matches(normalize(t),seen_taxa,n=4) for t in taxa}
            raise ValueError(f'No taxon matches with queried-rank confidence >=0.80 for {taxa}. Closest available names: {suggestions}. Only ranks present in SINTAX can be matched.')
        # Atomic exclusive publication; never replace an existing FASTA.
        os.link(temporary,output)
    print(f'Wrote {total:,} sequences to {output}',flush=True)
    return total


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--taxon',dest='taxa',action='append',required=True,help='Exact reported name, label, or taxid with queried-rank confidence >=0.80; repeat to select a union')
    parser.add_argument('--folders',nargs='+',required=True)
    parser.add_argument('--analysis-dir',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--merged-dir',type=Path)
    parser.add_argument('--sintax-root',type=Path)
    parser.add_argument('--suffix',default='sintaxSP')
    export_reads(**vars(parser.parse_args()))

if __name__=='__main__': main()
