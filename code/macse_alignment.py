"""Run MACSE in an isolated workspace and publish validated, traceable alignments."""
from pathlib import Path
import tempfile
import re
import subprocess


def align_macse(runner, retained, alignment, genetic_code=5, memory="2g"):
    """Align coding consensuses, trim ends, then restore tracking headers for IQ-TREE."""
    from postprocess_sintax import fasta, write_record, valid_alignment, tsv
    if not re.fullmatch(r'[1-9][0-9]*[mMgG]', memory):
        raise ValueError('MACSE memory must be a positive size such as 32g or 4096m')
    alignment = Path(alignment)
    prefix = str(alignment).removesuffix('.aln.fasta')
    def run_macse(command, log):
        # Surface useful diagnostic paths instead of only a long subprocess argument list.
        try:
            runner.run(command, log)
        except subprocess.CalledProcessError as error:
            runner.note(f'MACSE failed (exit {error.returncode}, Java heap limit {memory}). '
                        f'See {log} and {runner.log}. If memory is exhausted, resume '
                        'with a larger --macse-memory value.')
            raise
    runner.note(f'MACSE Java maximum heap: {memory}')
    originals = list(fasta(retained))
    labels = {f'seq{i:09d}': header for i, (header, _) in enumerate(originals)}
    # Temporary short IDs avoid Java/FASTA whitespace handling changing original headers.
    with tempfile.TemporaryDirectory(prefix='.macse_', dir=alignment.parent) as directory:
        work = Path(directory)
        source = work/'input.fasta'
        with source.open('w') as handle:
            for identifier, (_, sequence) in zip(labels, originals):
                write_record(handle, identifier, sequence)
        nt, aa, trimmed, report = (work/name for name in ('aligned.fasta', 'aligned.aa.fasta', 'trimmed.fasta', 'trimmed.csv'))
        run_macse(['macse', f'-Xmx{memory}', '-prog', 'alignSequences', '-gc_def', genetic_code,
                    '-gap_ext', '2.0', '-gap_op', '10', '-out_AA', aa, '-out_NT', nt,
                    '-seq', source], Path(prefix+'.macse_align.stdout.log'))
        # MACSE 2.07 accepts a fraction: 0.7 means 70%; literal 70 is rejected.
        run_macse(['macse', f'-Xmx{memory}', '-prog', 'trimAlignment', '-min_percent_NT_at_ends', '0.7',
                    '-out_NT', trimmed, '-out_trim_info', report, '-align', nt], Path(prefix+'.macse_trim.stdout.log'))
        # Validate every output before publishing; trimming may legitimately remove bases.
        restored = {}
        for path in (nt, aa, trimmed):
            records = list(fasta(path))
            ids = [header.split()[0] for header, _ in records]
            if len(ids) != len(labels) or set(ids) != set(labels):
                raise ValueError(f'MACSE changed the sequence set: {path.name}')
            if any(not sequence for _, sequence in records) or len({len(s) for _, s in records}) != 1:
                raise ValueError(f'MACSE produced an empty or unequal-length alignment: {path.name}')
            restored[path] = [(labels[h.split()[0]], s) for h, s in records]
        if not report.is_file() or not report.stat().st_size:
            raise ValueError('MACSE trim report is missing or empty')
        # Keep raw MACSE alignments, including frameshift symbols, for auditing.
        for source_path, suffix in ((nt, '.macse.aln.fasta'), (aa, '.aa.aln.fasta'), (trimmed, '.macse.trimmed.aln.fasta')):
            with Path(prefix+suffix).open('x') as handle:
                for header, sequence in restored[source_path]:
                    write_record(handle, header, sequence)
        replacements = []
        with alignment.open('x') as handle:
            for header, sequence in restored[trimmed]:
                # MACSE's ! denotes a frameshift; N makes this site missing data for IQ-TREE.
                sequence = sequence.upper()
                if set(sequence) - set('ACGTRYSWKMBDHVN-!'):
                    raise ValueError(f'Unexpected MACSE nucleotide symbols in {header}')
                replacements.append((header.split()[0], sequence.count('!')))
                write_record(handle, header, sequence.replace('!', 'N'))
        Path(prefix+'_trimmed.csv').write_bytes(report.read_bytes())
        tsv(Path(prefix+'.macse_ids.tsv'), ['macse_id', 'original_header'], labels.items())
        tsv(Path(prefix+'.frameshifts.tsv'), ['sequence_id', 'frameshifts_replaced_with_N'], replacements)
    valid_alignment(alignment, retained, allow_trimmed=True)
