"""Shared FASTA/FASTQ/TSV readers and writers used by both modules."""
import csv
import gzip
from pathlib import Path
import re

# fastplong names the pieces of a split chimeric read chimera-segment-N-<raw read ID>
# (older versions: split-by-adapter-left/right-<raw read ID>).
SPLIT_PREFIX = re.compile(r'^(?:split-by-adapter-(?:left|right)|chimera-segment-[0-9]+)-')


def raw_read_id(read_id):
    """ID of the raw ONT read a cleaned read (or a piece of it) comes from."""
    return SPLIT_PREFIX.sub('', read_id)


def open_text(path, mode='rt'):
    """Open plain or gzip-compressed text by file name."""
    path = Path(path)
    return gzip.open(path, mode) if path.suffix == '.gz' else path.open(mode.replace('t', ''))


def fasta(path):
    """Stream multiline FASTA records as (full header, sequence)."""
    header, parts = None, []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(parts)
                header, parts = line[1:], []
            else:
                if header is None:
                    raise ValueError(f'Invalid FASTA: {path}')
                parts.append(line)
    if header is not None:
        yield header, ''.join(parts)


def fastq_records(path):
    """Stream four-line FASTQ records as (read ID, sequence)."""
    with open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().rstrip('\r\n')
            plus = handle.readline()
            quality = handle.readline().rstrip('\r\n')
            if not header.startswith('@') or not plus.startswith('+') or len(sequence) != len(quality):
                raise ValueError(f'Invalid FASTQ record in {path}: {header[:100]!r}')
            yield header[1:].split()[0], sequence


def write_record(handle, header, sequence):
    """Write a FASTA record with consistent 80-base wrapping."""
    handle.write(f'>{header}\n')
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start:start + 80] + '\n')


def tsv(path, columns, rows):
    """Write a tab-delimited table with one header row (gzip-compressed if the name ends in .gz)."""
    with open_text(path, 'wt') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(columns)
        writer.writerows(rows)


def read_tsv(path):
    """Stream the rows of a tab-delimited table (plain or .gz) as dictionaries."""
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter='\t')
