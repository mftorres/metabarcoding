#!/usr/bin/env python3
"""Stream SINTAX TSVs into one provenance-preserving table; no classifier runs."""
import argparse
import csv
from pathlib import Path
import re

RANKS = dict(d='domain', k='kingdom', p='phylum', c='class', o='order',
             f='family', g='genus', s='species', t='strain')
BASE_COLUMNS = ['source_folder', 'barcode', 'source_name', 'source_file',
                'query_header', 'readid', 'raw_ont_read_id', 'strand', 'taxonomy',
                'taxlabel', 'taxonomy_lowest', 'taxonomy_lowestrank',
                'taxonomy_lowestboot', 'taxlabel_lowest', 'taxlabel_lowestrank',
                'sintax_taxid', 'sintax_name', 'classified']
ASSIGNMENT_COLUMNS = ['analysis_cutoff', 'assigned_rank_code', 'assigned_rank', 'assigned_label', 'assigned_name', 'assigned_taxid', 'assigned_confidence', 'taxonomy_status']
COLUMNS = BASE_COLUMNS + ASSIGNMENT_COLUMNS + [col for rank in RANKS.values()
                          for col in (rank, rank + '_bootstrap', 'accepted_' + rank, 'accepted_' + rank + '_confidence')]


def discover_outputs(inputs, suffix='sintaxSP'):
    """Directories select exact sample_suffix.tsv names, excluding derived subsets.

    Explicit TSV paths allow nonstandard names and different classifier runs.
    Full paths are retained to distinguish equal sample names from separate runs.
    """
    paths = set()
    for value in inputs:
        path = Path(value).resolve()
        if path.is_file():
            paths.add(path)
        elif path.is_dir():
            paths.update(p for p in path.rglob(f'*_{suffix}.tsv')
                         if p.name == f'{p.parent.name}_{suffix}.tsv')
        else:
            raise FileNotFoundError(path)
    if not paths:
        raise FileNotFoundError(f'No *_{suffix}.tsv classifier outputs found in {inputs}')
    return sorted(paths)


def parse_taxonomy(value):
    result = {}
    for token in value.strip().split(','):
        if not token:
            continue
        match = re.fullmatch(r'([dkpcofgst]):(.+?)(?:\(([0-9.]+)\))?', token.strip())
        if not match:
            raise ValueError(f'Malformed taxonomy token: {token!r}')
        rank, name, confidence = match.groups()
        if rank in result:
            raise ValueError(f'Duplicate rank: {rank}')
        boot = float(confidence) if confidence is not None else None
        if boot is not None and not 0 <= boot <= 1:
            raise ValueError(f'Bootstrap outside [0, 1]: {confidence}')
        result[rank] = (name, boot)
    return result


def taxon_parts(value):
    """Only interpret a trailing numeric MIDORI suffix as an NCBI taxid."""
    match = re.fullmatch(r'(.+)_([0-9]+)', value)
    return (match[1], match[2]) if match else (value, '')


def parse_output(path, cutoff=None):
    if cutoff is not None and not 0 <= cutoff <= 1:
        raise ValueError("cutoff must be in [0, 1]")
    path = Path(path).resolve()
    with path.open(encoding='utf-8', newline='') as handle:
        for line_number, fields in enumerate(csv.reader(handle, delimiter='\t'), 1):
            if not fields:
                continue
            if len(fields) not in (3, 4) or not fields[0].strip():
                raise ValueError(f'{path}:{line_number}: expected 3 or 4 SINTAX columns')
            header, taxonomy, strand = fields[:3]
            if strand not in ('', '+', '-', '*'):
                raise ValueError(f'{path}:{line_number}: invalid strand {strand!r}')
            # Preserve column four; an explicit analysis cutoff recomputes accepted ranks from column two.
            label = fields[3] if len(fields) == 4 else ''
            try:
                candidates, accepted = parse_taxonomy(taxonomy), parse_taxonomy(label)
            except ValueError as error:
                raise ValueError(f'{path}:{line_number}: {error}') from error
            if cutoff is not None:
                accepted = {code: (name, boot) for code, (name, boot) in candidates.items()
                            if boot is not None and boot >= cutoff}
            row = dict.fromkeys(COLUMNS, '')
            row['analysis_cutoff'] = '' if cutoff is None else cutoff
            row['taxonomy_status'] = 'taxafail' 
            readid = header.split()[0]
            raw_id = re.sub(r'^split-by-adapter-(?:left|right)-', '', readid)
            parent = re.search(r'(?:^|\s)parent_read_id=([^\s]+)', header)
            if parent:
                raw_id = parent[1]
            row.update(source_folder=path.parent.name, barcode=path.parent.name,
                       source_name=path.name, source_file=str(path), query_header=header,
                       readid=readid, raw_ont_read_id=raw_id, strand=strand,
                       taxonomy=taxonomy, taxlabel=label, classified=bool(accepted))
            for code, rank in RANKS.items():
                name, boot = candidates.get(code, ('', None))
                row[rank] = name
                row[rank + '_bootstrap'] = '' if boot is None else boot
                row['accepted_' + rank] = accepted.get(code, ('', None))[0]
                if code in accepted:
                    row['accepted_' + rank + '_confidence'] = candidates.get(code, ('', ''))[1]
            if candidates:
                code = next(code for code in reversed(RANKS) if code in candidates)
                name, boot = candidates[code]
                row.update(taxonomy_lowest=name, taxonomy_lowestrank=code,
                           taxonomy_lowestboot='' if boot is None else boot)
            if accepted:
                code = next(code for code in reversed(RANKS) if code in accepted)
                name = accepted[code][0]
                clean_name, taxid = taxon_parts(name)
                row.update(taxlabel_lowest=name, taxlabel_lowestrank=code,
                           sintax_taxid=taxid, sintax_name=clean_name)
                row.update(assigned_rank_code=code, assigned_rank=RANKS[code], assigned_label=name,
                           assigned_name=clean_name, assigned_taxid=taxid,
                           assigned_confidence=candidates.get(code, ('', ''))[1],
                           taxonomy_status='taxapass' if code in ('g', 's') else 'taxafail')
            yield row


def load_outputs(inputs, suffix='sintaxSP', cutoff=None):
    import pandas as pd
    rows = (row for path in discover_outputs(inputs, suffix) for row in parse_output(path, cutoff=cutoff))
    frame = pd.DataFrame(rows, columns=COLUMNS)
    for col in ['taxonomy_lowestboot'] + [r + '_bootstrap' for r in RANKS.values()]:
        frame[col] = pd.to_numeric(frame[col], errors='coerce')
    frame['sintax_taxid'] = pd.to_numeric(frame['sintax_taxid'], errors='coerce').astype('Int64')
    return frame


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs', nargs='+', type=Path)
    parser.add_argument('--suffix', default='sintaxSP')
    parser.add_argument('--cutoff', type=float, default=0.8)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    paths = discover_outputs(args.inputs, args.suffix)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive creation prevents accidental replacement, including an input TSV.
    with args.output.open('x', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, COLUMNS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        count = 0
        for path in paths:
            for row in parse_output(path, cutoff=args.cutoff):
                writer.writerow(row)
                count += 1
    print(f'Parsed {count:,} records from {len(paths)} files into {args.output}')


if __name__ == '__main__':
    main()
