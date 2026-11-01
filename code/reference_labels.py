"""Stable reference IDs with tool-compatible taxon suffixes."""
import csv
from pathlib import Path
import re
from genbank_xml import iter_genbank_records

TAXON_SEPARATOR = '__taxon__'


def reference_id(identifier):
    """Recover the source identifier from either old or annotated headers."""
    return identifier.split(TAXON_SEPARATOR, 1)[0]


def labelled_reference_id(identifier, taxon):
    # Avoid whitespace and punctuation interpreted by Newick or alignment tools.
    label = re.sub(r'[^A-Za-z0-9_.-]+', '_', taxon.strip()).strip('_') or 'unassigned'
    return reference_id(identifier) + TAXON_SEPARATOR + label


def reference_names(path):
    names = {}
    if Path(path).suffix == '.tsv':
        with Path(path).open() as handle:
            return {row['sequence_id']: row['taxon'] for row in csv.DictReader(handle, delimiter='\t')}
    for record in iter_genbank_records(path):
        organism = record.findtext('GBSeq_organism')
        if organism:
            for field in ('GBSeq_accession-version', 'GBSeq_primary-accession'):
                identifier = record.findtext(field)
                if identifier:
                    names[identifier] = organism
    return names
