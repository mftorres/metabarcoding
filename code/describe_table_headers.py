#!/usr/bin/env python3
"""Validate analysis TSV headers and document their columns in TABLE_HEADERS.md."""
import argparse
import csv
from collections import defaultdict
from pathlib import Path


def describe_headers(directory):
    directory=Path(directory)
    groups=defaultdict(list)
    paths=sorted(directory.rglob('*.tsv'))
    if not paths: raise ValueError(f'No TSV tables in {directory}')
    for path in paths:
        with path.open(encoding='utf-8',newline='') as handle:
            reader=csv.reader(handle,delimiter='\t')
            columns=next(reader,None)
            if not columns or any(not column for column in columns) or len(set(columns))!=len(columns):
                raise ValueError(f'Missing, blank or duplicate column headers: {path}')
            if 'source_folder' not in columns and 'rank' not in columns:
                raise ValueError(f'Unexpected table schema: {path}')
            row=next(reader,None)
            if row is not None and len(row)!=len(columns):
                raise ValueError(f'Header/data column mismatch: {path}')
        if path.name=='species_by_folder.tsv':
            key=('species_by_folder.tsv',tuple(columns[:2]+columns[-2:]))
            groups[key].append((path,len(columns)-4))
        else:
            groups[('tables',tuple(columns))].append((path,None))
    text=['# Output table headers\n',
          'Every `.tsv` is tab-delimited, with one header row as the first line. Import it with `sep="\\t", header=0`.\n',
          'Primer columns contain mapped PCR names: dictionary keys `forward_tagged_primer_1` / `reverse_tagged_primer_1` map to values `bc1002` / `bc1100rc`. Original-name columns are not exported.\n',
          '`forward_PCR_primer` and `reverse_PCR_primer` are read/primer-pair fields. `forward_primers_detected` and `reverse_primers_detected` are unique semicolon-separated lists.\n']
    for (kind,columns),members in groups.items():
        examples=', '.join(f'`{str(p.relative_to(directory))}`' for p,_ in members[:2])
        if len(members)>2: examples+=f' (same headers in {len(members)} tables)'
        text+=['\n## '+examples+'\n','| Position | Column |\n|---:|---|\n']
        if kind=='species_by_folder.tsv':
            count=members[0][1]
            text += [f'| 1 | `{columns[0]}` |\n',f'| 2 | `{columns[1]}` |\n',
                     f'| 3–{count+2} | {count:,} dynamic species-label columns, each containing confident read counts |\n',
                     f'| {count+3} | `{columns[2]}` |\n',f'| {count+4} | `{columns[3]}` |\n']
        else:
            text += [f'| {i} | `{col}` |\n' for i,col in enumerate(columns,1)]
    target=directory/'TABLE_HEADERS.md'
    target.write_text(''.join(text),encoding='utf-8')
    print(f'Validated header rows in {len(paths)} TSV tables; wrote {target}',flush=True)
    return len(paths)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    describe_headers(parser.parse_args().directory)

if __name__=='__main__': main()
