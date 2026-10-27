#!/usr/bin/env python3
"""Inventory/remove obsolete consensus, clustering and previous taxonomy outputs.

Default is a TSV plan. --apply deletes only inventoried files after verifying their
size and modification time. Raw/cleaned reads, primer results, reference databases,
and read-level SINTAX results are outside the removal rules.
"""
import argparse
import csv
from datetime import datetime,timezone
from pathlib import Path


def reason_for(path,root):
    rel=path.relative_to(root)
    parts=rel.parts
    name=path.name.lower()
    if not parts: return None
    if parts[0]=='data':
        if len(parts)>1 and parts[1]=='kraken_parsed': return 'previous classifier parsed output'
        if 'consensus' in parts and name.endswith(('.fasta','.fa','.fas','.fna','.fasta.gz','.fa.gz')):
            return 'consensus FASTA'
        if 'nonumt' in name or name=='no_consensus.done': return 'retired nucleotide-translation/noNUMT output'
        if name.endswith('.kraken') or name.endswith('.fasta.report'): return 'previous classifier output'
        if 'consensus' in name and 'sintax' in name and path.suffix=='.tsv': return 'classification of retired consensus sequences'
        if len(parts)>1 and parts[1].startswith('mergedseq'):
            if 'filtered' in parts:
                index=parts.index('filtered')
                if index+1<len(parts) and parts[index+1]=='combined': return 'retired clustering output'
                if 'amplicon_sorter' in name or 'kraken2' in name or name in ('consensusfile.fasta','results.csv'):
                    return 'retired clustering/classifier output or marker'
            if any(token in name for token in ('amplicon_sorter','read_primer_group_taxid','target_species_detection',
                'target_species_taxids_below_50pct','taxids_across_50pct','top10_taxa_per_forward_tag',
                'unsorted_taxa','expected_species_not_detected','primer_detection_classification_summary')):
                return 'superseded group-based taxonomy table'
            if 'figures' in parts and any(token in name for token in ('target_detection','target_species_taxids',
                'top10_taxids','unique_taxids','insecta_reads_vs_taxid','reads_per_forward_primer',
                'reads_per_reverse_primer','primer_detection_classification')):
                return 'superseded group-based taxonomy figure'
    if len(parts)>1 and parts[:2]==('code','executed_notebooks') and len(parts)==3:
        if name.startswith(('kraken_parser','evaluate_tagtracking_taxid_reduced','sankey_target_species_detection')):
            return 'superseded executed taxonomy notebook'
    if len(parts)>1 and parts[:2]==('code','.ipynb_checkpoints') and name.startswith('evaluate_tagtracking'):
        return 'superseded notebook checkpoint'
    return None


def inventory(root):
    for base in ('data','code/executed_notebooks','code/.ipynb_checkpoints'):
        directory=root/base
        if not directory.exists(): continue
        for path in directory.rglob('*'):
            if not path.is_file() or path.is_symlink(): continue
            reason=reason_for(path,root)
            if reason:
                stat=path.stat()
                yield dict(path=str(path.relative_to(root)),bytes=stat.st_size,mtime_ns=stat.st_mtime_ns,reason=reason)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--apply',action='store_true')
    args=parser.parse_args();root=args.root.resolve()
    if not args.apply:
        rows=list(inventory(root));args.manifest.parent.mkdir(parents=True,exist_ok=True)
        with args.manifest.open('x',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=['path','bytes','mtime_ns','reason'],delimiter='\t',lineterminator='\n')
            writer.writeheader();writer.writerows(rows)
        print(f'Planned {len(rows):,} files; {sum(r["bytes"] for r in rows)/1e9:.3f} GB. Manifest: {args.manifest}')
        return
    with args.manifest.open() as handle: rows=list(csv.DictReader(handle,delimiter='\t'))
    # Validate the full deletion list before removing any files.
    for row in rows:
        path=root/row['path']
        if not path.resolve().is_relative_to(root) or not reason_for(path,root) or path.is_symlink():
            raise ValueError(f'Unsafe or unrecognised path: {path}')
        stat=path.stat()
        if stat.st_size!=int(row['bytes']) or stat.st_mtime_ns!=int(row['mtime_ns']):
            raise ValueError(f'File changed after inventory: {path}')
    parents=set()
    for row in rows:
        path=root/row['path'];parents.add(path.parent);path.unlink()
    for parent in sorted(parents,key=lambda p:len(p.parts),reverse=True):
        try: parent.rmdir()
        except OSError: pass
    log=args.manifest.with_name(args.manifest.stem+'_completed.tsv')
    with log.open('x',newline='') as handle:
        writer=csv.writer(handle,delimiter='\t',lineterminator='\n')
        writer.writerow(['completed_utc','deleted_files','deleted_bytes'])
        writer.writerow([datetime.now(timezone.utc).isoformat(),len(rows),sum(int(r['bytes']) for r in rows)])
    print(f'Removed {len(rows):,} inventoried obsolete files. Record: {log}')

if __name__=='__main__':main()
