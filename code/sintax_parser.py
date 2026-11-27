"""Read VSEARCH SINTAX output: one line per classified sequence.

Each line has the sequence name, the taxonomy with a confidence per rank, for example
k:Eukaryota_2759(1.00),p:Arthropoda_6656(1.00),...,s:Phaonia mystica_2882605(0.93),
the strand, and SINTAX's own cutoff column (not used: the cutoff is applied later).
"""
import re

# MIDORI labels these ranks; the trailing number of a name is its NCBI taxid.
RANKS = dict(k='kingdom', p='phylum', c='class', o='order', f='family', g='genus', s='species')


def taxon_parts(label):
    """Split a MIDORI label such as 'Phaonia mystica_2882605' into name and NCBI taxid."""
    match = re.fullmatch(r'(.+)_([0-9]+)', label)
    if match:
        return match[1], match[2]
    return label, ''


def parse_taxonomy(value):
    """{rank name: (label, confidence)} from a SINTAX taxonomy string."""
    result = {}
    for token in value.split(','):
        if not token:
            continue
        match = re.fullmatch(r'([a-z]):(.+)\(([0-9.]+)\)', token.strip())
        if not match or match[1] not in RANKS:
            raise ValueError(f'Unexpected taxonomy item {token!r}; the pipeline expects the MIDORI database')
        result[RANKS[match[1]]] = (match[2], float(match[3]))
    return result


def parse_output(path):
    """Yield (sequence name, strand, taxonomy string, {rank: (label, confidence)}) per line."""
    with open(path, encoding='utf-8') as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 3 or not fields[0]:
                raise ValueError(f'{path}, line {number}: expected at least 3 tab-separated SINTAX columns')
            name, taxonomy, strand = fields[0].split()[0], fields[1], fields[2]
            try:
                ranks = parse_taxonomy(taxonomy)
            except ValueError as error:
                raise ValueError(f'{path}, line {number}: {error}') from error
            yield name, strand, taxonomy, ranks
