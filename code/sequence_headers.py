"""FASTA headers for local reads exported to Postprocessing.

Format: folder-forwardtag-reversetag-trustedtaxon-contamination-readid. Fields are
percent-encoded so they never contain '-'; a missing tag is written as 'none'.
"""
import re
from urllib.parse import quote, unquote

FLAGS = ('contamination', 'nocontamination')


def header_component(value):
    return quote(re.sub(r'\s+', '_', value), safe='_.').replace('-', '%2D')


def export_header(folder, forward, reverse, assigned_label, contamination, readid):
    return '-'.join([header_component(folder), header_component(forward or 'none'),
                     header_component(reverse or 'none'), header_component(assigned_label or 'unassigned'),
                     contamination, readid])


def local_fields(identifier):
    """Split a local header into its six fields, or return None for reference sequences."""
    parts = identifier.split()[0].split('-', 5)
    return parts if len(parts) == 6 and parts[4] in FLAGS else None


def local_annotations(header):
    """(folder, forward tag, reverse tag, trusted taxon) of an exported local read."""
    parts = local_fields(header)
    if parts is None:
        raise ValueError(f'Expected an annotated local read header: {header}')
    return tuple(unquote(value) for value in parts[:4])
