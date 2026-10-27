"""Taxonomy lookup from MIDORI FASTA headers for notebook annotation."""
from functools import lru_cache
from sintax_parser import parse_taxonomy, taxon_parts, RANKS

@lru_cache(maxsize=1)
def catalog(path):
    names, parents, ranks = {}, {}, {}
    with open(path, encoding='utf-8') as handle:
        for line in handle:
            if not line.startswith('>') or 'tax=' not in line:
                continue
            taxonomy = parse_taxonomy(line.split('tax=', 1)[1].strip().rstrip(';'))
            previous = None
            for code in RANKS:
                if code not in taxonomy:
                    continue
                name, value = taxon_parts(taxonomy[code][0])
                if not value:
                    continue
                taxid = int(value)
                # Some SP entries repeat a species taxid at synthetic genus ranks.
                # Prefer the deepest rank and never create a self-parent cycle.
                names[taxid] = name
                ranks[taxid] = RANKS[code]
                if previous != taxid:
                    parents[taxid] = previous if previous is not None else taxid
                previous = taxid
    return names, parents, ranks

def read_scientific_names(reference, wanted_taxids):
    names, _, _ = catalog(str(reference))
    return {taxid: names[taxid] for taxid in wanted_taxids if taxid in names}

def scientific_name_taxids(reference, wanted_name_keys):
    import re
    names, _, _ = catalog(str(reference))
    matches = {name: [] for name in wanted_name_keys}
    for taxid, name in names.items():
        key = re.sub(r'[\s_]+', '_', name.strip().lower()).strip('_')
        if key in matches:
            matches[key].append(taxid)
    return matches

def lineage_subset(reference, start_taxids):
    _, parents, ranks = catalog(str(reference))
    return parents, ranks
