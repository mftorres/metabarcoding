"""Parse VSEARCH SINTAX tab-delimited output into one row per classified read."""
import csv
from pathlib import Path
import re

RANKS = dict(d='domain', k='kingdom', p='phylum', c='class', o='order',
             f='family', g='genus', s='species', t='strain')
BASE_COLUMNS = ['source_folder', 'query_header', 'readid', 'raw_ont_read_id', 'strand', 'taxonomy',
                'taxlabel', 'taxonomy_lowest', 'taxonomy_lowestrank', 'taxonomy_lowestboot', 'classified']
ASSIGNMENT_COLUMNS = ['analysis_cutoff', 'assigned_rank_code', 'assigned_rank', 'assigned_label', 'assigned_name',
                      'assigned_taxid', 'assigned_confidence', 'taxonomy_status']
COLUMNS = BASE_COLUMNS + ASSIGNMENT_COLUMNS + [col for rank in RANKS.values()
                          for col in (rank, rank + '_bootstrap', 'accepted_' + rank, 'accepted_' + rank + '_confidence')]
SPLIT_PREFIX = re.compile(r'^(?:split-by-adapter-(?:left|right)|chimera-segment-[0-9]+)-')


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


def raw_read_id(readid):
    """ID of the raw read. fastplong names pieces of split chimeric reads
    chimera-segment-N-<raw ID> (older versions: split-by-adapter-left/right-<raw ID>)."""
    return SPLIT_PREFIX.sub('', readid)


def parse_output(path, cutoff):
    """Yield one row per SINTAX query; accepted ranks are those with confidence >= cutoff."""
    if not 0 <= cutoff <= 1:
        raise ValueError('cutoff must be in [0, 1]')
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
            # Column four (SINTAX's own cutoff) is preserved; acceptance is recomputed here.
            label = fields[3] if len(fields) == 4 else ''
            try:
                candidates = parse_taxonomy(taxonomy)
                parse_taxonomy(label)
            except ValueError as error:
                raise ValueError(f'{path}:{line_number}: {error}') from error
            accepted = {code: (name, boot) for code, (name, boot) in candidates.items()
                        if boot is not None and boot >= cutoff}
            row = dict.fromkeys(COLUMNS, '')
            readid = header.split()[0]
            row.update(source_folder=path.parent.name, query_header=header, readid=readid,
                       raw_ont_read_id=raw_read_id(readid), strand=strand, taxonomy=taxonomy,
                       taxlabel=label, classified=bool(accepted), analysis_cutoff=cutoff,
                       taxonomy_status='taxafail')
            for code, rank in RANKS.items():
                name, boot = candidates.get(code, ('', None))
                row[rank] = name
                row[rank + '_bootstrap'] = '' if boot is None else boot
                if code in accepted:
                    row['accepted_' + rank] = name
                    row['accepted_' + rank + '_confidence'] = boot
            if candidates:
                code = next(code for code in reversed(RANKS) if code in candidates)
                name, boot = candidates[code]
                row.update(taxonomy_lowest=name, taxonomy_lowestrank=code,
                           taxonomy_lowestboot='' if boot is None else boot)
            if accepted:
                code = next(code for code in reversed(RANKS) if code in accepted)
                name, boot = accepted[code]
                clean_name, taxid = taxon_parts(name)
                row.update(assigned_rank_code=code, assigned_rank=RANKS[code], assigned_label=name,
                           assigned_name=clean_name, assigned_taxid=taxid, assigned_confidence=boot,
                           taxonomy_status='taxapass' if code in ('g', 's') else 'taxafail')
            yield row
