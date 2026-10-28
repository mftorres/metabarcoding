"""Isolate kc-align fixed filenames and restore consensus tracking headers."""
from pathlib import Path
import tempfile
import subprocess
import sys


def validate_reference(path, table=5):
    """Require one unambiguous, in-frame coding reference before starting a run."""
    from postprocess_sintax import fasta
    from Bio.Seq import Seq
    records = list(fasta(path))
    if len(records) != 1:
        raise ValueError('--kcalign-reference must contain exactly one nucleotide sequence')
    header, sequence = records[0]
    sequence = sequence.upper()
    if not sequence or len(sequence) % 3 or set(sequence) - set('ACGT'):
        raise ValueError('Trusted reference must be unambiguous A/C/G/T, ungapped and in-frame (length divisible by 3)')
    if '*' in str(Seq(sequence).translate(table=table))[:-1]:
        raise ValueError('Trusted reference has internal stop codons under the chosen translation table')
    return header, sequence


def align_kcalign(runner, retained, alignment, reference, genetic_code=5, threads=6):
    """Run mixed mode on selected consensuses; audit reads rejected by the aligner."""
    from postprocess_sintax import fasta, write_record, tsv, validate_outputs
    _, refseq = validate_reference(reference, genetic_code)
    alignment = Path(alignment)
    prefix = str(alignment).removesuffix('.aln.fasta')
    records = list(fasta(retained))
    labels = {f'query{i:09d}': header for i, (header, _) in enumerate(records)}
    # kc-align adds its reference to the alignment; remove only this synthetic guide ID.
    guide = 'trusted_alignment_guide'
    with tempfile.TemporaryDirectory(prefix='.kcalign_', dir=alignment.parent) as directory:
        work = Path(directory).resolve()
        with (work/'input.fasta').open('w') as handle:
            for identifier, (_, sequence) in zip(labels, records):
                write_record(handle, identifier, sequence)
        with (work/'reference.fasta').open('w') as handle:
            write_record(handle, guide, refseq)
        command = [sys.executable, Path(__file__).with_name('kcalign_cli.py').resolve(),
                   '--mode', 'mixed', '--reference', work/'reference.fasta',
                   '--sequences', work/'input.fasta', '--table', genetic_code,
                   '--keep', '--threads', threads]
        # Fixed tool filenames are confined to this directory; original data stay unchanged.
        runner.run(command, Path(prefix+'.kcalign.stdout.log'), cwd=work)
        outputs = {}
        for name in ('codon_align.fasta', 'pre_align.fasta'):
            validate_outputs([work/name])
            rows = list(fasta(work/name))
            unknown = {h.split()[0] for h, _ in rows} - set(labels) - {guide}
            if unknown:
                raise ValueError(f'Unexpected kc-align IDs: {unknown}')
            outputs[name] = [(labels[h.split()[0]], s) for h, s in rows if h.split()[0] != guide]
        aligned = outputs['codon_align.fasta']
        if not aligned or len({len(s) for _, s in aligned}) != 1:
            raise ValueError('kc-align produced an empty or unequal-length alignment')
        if any(set(s.upper()) - set('ACGTRYSWKMBDHVN-') for _, s in aligned):
            raise ValueError('kc-align produced invalid nucleotide alignment characters')
        kept = {h for h, _ in aligned}
        # Mixed mode can reject reads; record every input outcome rather than inventing tips.
        tsv(Path(prefix+'.kcalign_selection.tsv'), ['sequence_id', 'status'],
            [(h.split()[0], 'aligned' if h in kept else 'rejected_by_kcalign') for h, _ in records])
        tsv(Path(prefix+'.kcalign_ids.tsv'), ['kcalign_id', 'original_header'], labels.items())
        if len(aligned) < 4:
            raise ValueError('Fewer than four consensus sequences survived kc-align; inspect its log and selection table')
        with Path(prefix+'.pre_align.aa.fasta').open('x') as handle:
            for h, s in outputs['pre_align.fasta']:
                write_record(handle, h, s)
        with alignment.open('x') as handle:
            for h, s in aligned:
                write_record(handle, h, s)
