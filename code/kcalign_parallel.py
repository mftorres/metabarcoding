"""Safe parallel mixed-mode screening using kc-align's serial per-read checks."""
from concurrent.futures import ProcessPoolExecutor
import os
from pathlib import Path
import tempfile
import kcalign
from Bio import SeqIO

# Keep the installed implementation: its serial rejection path preserves read IDs.
_original_create_lists = kcalign.create_lists


def screen_read(task):
    """Isolate fixed temporary filenames and run the original checks on one read."""
    record, reference_protein, reference_nt, table, distance = task
    original_directory = Path.cwd()
    with tempfile.TemporaryDirectory(prefix='kc_read_') as directory:
        try:
            os.chdir(directory)
            SeqIO.write([record], 'read.fasta', 'fasta')
            return _original_create_lists('read.fasta', reference_protein,
                                          dict(reference_nt), 0, table, distance, 1)
        finally:
            os.chdir(original_directory)


def merge_results(results, reference_nt):
    """Combine accepted records and rejected IDs without interpreting rejection as data."""
    proteins, descriptions, identifiers, rejected = [], [], [], []
    nucleotides = dict(reference_nt)
    for seqs, names, ids, sequences, errors in results:
        proteins.extend(seqs)
        descriptions.extend(names)
        identifiers.extend(ids)
        rejected.extend(errors)
        for identifier in ids:
            nucleotides[identifier] = sequences[identifier]
    return proteins, descriptions, identifiers, nucleotides, rejected


def create_lists(reads, reference_protein, reference_nt, join, table, distance, cores):
    """Replace only the broken mixed-mode parallel dispatcher, not biological filters."""
    if join != 0 or int(cores) == 1:
        return _original_create_lists(reads, reference_protein, reference_nt, join, table, distance, cores)
    records = [record for record in SeqIO.parse(reads, 'fasta') if record.id not in reference_nt]
    tasks = [(record, reference_protein, reference_nt, table, distance) for record in records]
    # Preserve the installed tool's three-thread budget per concurrently screened read.
    # A single pool also avoids launching a new process pool for every pair of reads.
    with ProcessPoolExecutor(max_workers=max(1, int(cores)//3)) as pool:
        return merge_results(pool.map(screen_read, tasks), reference_nt)
