#!/usr/bin/env python3
"""Validate analysis TSV headers and document their columns in TABLE_HEADERS.md."""
import argparse
import csv
from collections import defaultdict
from pathlib import Path


def describe_headers(directory):
    directory = Path(directory)
    groups = defaultdict(list)
    paths = sorted(directory.rglob('*.tsv'))
    if not paths:
        raise ValueError(f'No TSV tables in {directory}')
    for path in paths:
        with path.open(encoding='utf-8', newline='') as handle:
            reader = csv.reader(handle, delimiter='\t')
            columns = next(reader, None)
            if not columns or any(not column for column in columns) or len(set(columns)) != len(columns):
                raise ValueError(f'Missing, blank or duplicate column headers: {path}')
            row = next(reader, None)
            if row is not None and len(row) != len(columns):
                raise ValueError(f'Header/data column mismatch: {path}')
        if path.name == 'species_by_folder.tsv':
            # Species columns differ per dataset; describe them as one dynamic block.
            groups[('species_by_folder', tuple(columns[:1] + columns[-2:]))].append((path, len(columns) - 3))
        else:
            groups[('table', tuple(columns))].append((path, None))
    text = ['# Output table headers\n\n',
            'Every `.tsv` is tab-delimited with one header row. In Python: `pandas.read_csv(path, sep="\\t")`.\n\n',
            '`forward_PCR_primer` / `reverse_PCR_primer` hold the tag name (for example `bc1002`, `bc1100rc`) '
            'when the tag was assigned; `forward_tag_status` / `reverse_tag_status` say why it is blank otherwise. '
            '`*_primers_detected` columns are unique semicolon-separated lists.\n']
    for (kind, columns), members in groups.items():
        examples = ', '.join(f'`{p.relative_to(directory)}`' for p, _ in members[:2])
        if len(members) > 2:
            examples += f' (same headers in {len(members)} tables)'
        text += ['\n## ' + examples + '\n\n', '| Position | Column |\n|---:|---|\n']
        if kind == 'species_by_folder':
            count = members[0][1]
            text += [f'| 1 | `{columns[0]}` |\n',
                     f'| 2–{count + 1} | {count:,} species columns, each holding confident read counts |\n',
                     f'| {count + 2} | `{columns[1]}` |\n', f'| {count + 3} | `{columns[2]}` |\n']
        else:
            text += [f'| {i} | `{col}` |\n' for i, col in enumerate(columns, 1)]
    target = directory/'TABLE_HEADERS.md'
    target.write_text(''.join(text), encoding='utf-8')
    print(f'Validated header rows in {len(paths)} TSV tables; wrote {target}', flush=True)
    return len(paths)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    describe_headers(parser.parse_args().directory)


if __name__ == '__main__':
    main()
