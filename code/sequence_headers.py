"""Shared versioned headers for annotated local nucleotide sequence exports."""
import csv
import json
from pathlib import Path
import re
from urllib.parse import quote

HEADER_SCHEMA = 'folder-forward-reverse-trustedtaxon-contamination-readid-v2'


def header_component(value):
    return quote(re.sub(r'\s+', '_', value), safe='_.').replace('-', '%2D')


def export_header(folder, forward, reverse, assigned_label, contamination, readid):
    return '-'.join([header_component(folder), header_component(forward), header_component(reverse),
                     header_component(assigned_label or 'unassigned'), contamination, readid])


def local_header_parts(identifier):
    parts = identifier.split('-', 5)
    if len(parts) == 6 and parts[4] in ('contamination', 'nocontamination'):
        return parts[0], parts[5]
    # Existing figures can still read the original five-field local identifiers.
    old = identifier.split('-', 4)
    if len(old) == 5 and old[2] in ('taxapass', 'taxafail') and old[3] in ('contamination', 'nocontamination'):
        return old[0], old[4]
    return None


def primer_annotations(analysis_dir, folder, queries):
    """Map cleaned query IDs to canonical forward/reverse tags using raw parents.

    queries maps cleaned query IDs to parent raw IDs. Prefer the completed read
    table; fall back to the saved best-hit table if that table is not available.
    """
    if not queries:
        return {}
    analysis_dir = Path(analysis_dir)
    params = json.loads((analysis_dir / 'analysis_parameters.json').read_text())
    mapping = json.loads((analysis_dir / 'primer_name_mapping.json').read_text())
    canonical = {tag for value in mapping.values() for tag in value.split(';')}
    def mapped(value):
        value = mapping.get(value, value)
        if value and not set(value.split(';')).issubset(canonical):
            raise ValueError(f'Unknown primer name in analysis: {value}')
        return value or 'none'
    table = analysis_dir / folder / f'{folder}_reads.tsv'
    found = {}
    if table.is_file():
        with table.open() as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if {'forward_PCR_primer', 'reverse_PCR_primer', 'readid'}.issubset(reader.fieldnames or []):
                for row in reader:
                    query = row['readid']
                    if query in queries:
                        if row.get('raw_ont_read_id', queries[query]) != queries[query]:
                            raise ValueError(f'Analysis parent-ID mismatch for {query}')
                        missing = row.get('primer_tracking_record_found', '').lower() == 'false'
                        found[query] = ('unknown', 'unknown') if missing else (
                            mapped(row['forward_PCR_primer']), mapped(row['reverse_PCR_primer']))
                        if len(found) == len(queries):
                            break
    unresolved = {parent for query, parent in queries.items() if query not in found}
    if unresolved:
        parents = {}
        dataset = Path(params['merged_dir']).name
        with Path(params['primer_hits']).open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                if row['merged_folder'] == dataset and row['barcode_folder'] == folder and row['seqID'] in unresolved:
                    if row['seqID'] in parents:
                        raise ValueError(f'Duplicate best-hit row for {folder}/{row["seqID"]}')
                    parents[row['seqID']] = (mapped(row['forward_patternName']), mapped(row['reverse_patternName']))
        for query, parent in queries.items():
            if query not in found:
                found[query] = parents.get(parent, ('unknown', 'unknown'))
    return found
