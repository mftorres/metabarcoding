#!/usr/bin/env python3
"""Module 2, Postprocessing: export the clusters of one taxon, add references, align and build a draft tree."""
import argparse
from collections import Counter, OrderedDict
import csv
from datetime import date
import hashlib
import json
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

from bold_references import database_files, extract_bold
from genbank_xml import iter_genbank_records
from reference_labels import TAXON_SEPARATOR, reference_names, labelled_reference_id, reference_id
from sample_folders import analysis_folders
from export_clusters import header_fields
from seqio import fasta, write_record, tsv

# Settings that name the run's files; they cannot change on --resume.
FROZEN = ('taxon', 'folders', 'analysis_dir', 'date', 'reference_source')
CODE = Path(__file__).resolve().parent


def select_main_taxon(source, summary, output):
    """Within each specimen (barcode + forward tag), keep the clusters of the identification
    with most reads (all of them when tied). Return the number of clusters kept."""
    reads, clusters = Counter(), Counter()
    for header, _ in fasta(source):
        folder, forward, taxon, _, read_count, _ = header_fields(header)
        reads[folder, forward, taxon] += read_count
        clusters[folder, forward, taxon] += 1
    largest = {}
    for (folder, forward, taxon), count in reads.items():
        largest[folder, forward] = max(largest.get((folder, forward), 0), count)
    kept = {key for key, count in reads.items() if count == largest[key[:2]]}
    tsv(summary, ['source_folder', 'forward_tag', 'taxon', 'cluster_count', 'read_count', 'kept_for_tree'],
        [(*key, clusters[key], reads[key], str(key in kept).lower()) for key in sorted(reads)])
    count = 0
    with Path(output).open('x') as handle:
        for header, sequence in fasta(source):
            if header_fields(header)[:3] in kept:
                write_record(handle, header, sequence)
                count += 1
    return count


# Exclude genome descriptions, not mitochondrial origin: COI itself is mitochondrial.
GENOME_TITLE_QUERY = ('"complete genome"[Title] OR "partial genome"[Title] OR '
                      'mitogenome[Title] OR mitogenomes[Title] OR '
                      '"mitochondrial genome"[Title] OR "mitochondrial genomes"[Title] OR '
                      '((mitochondrion[Title] OR mitochondria[Title]) AND "complete sequence"[Title])')


def reference_query(taxid, include_mitogenomes=False):
    """Build the COI Entrez query with the chosen whole-genome exclusion policy."""
    query = (f'txid{taxid}[Organism:exp] AND (COI[Gene Name] OR CO1[Gene Name] OR '
             'COX1[Gene Name] OR "COI"[Title] OR "cytochrome oxidase i"[Title] OR '
             '"cytochrome oxidase 1"[Title] OR "cytochrome c oxidase subunit I"[Title] OR '
             '"cytochrome c oxidase subunit 1"[Title])')
    if not include_mitogenomes:
        query += f' NOT ({GENOME_TITLE_QUERY})'
    return query


def genome_description(header):
    """Recognise whole-genome titles without excluding ordinary mitochondrial genes."""
    title = header.partition(' ')[2].lower()
    return bool(re.search(r'\b(?:complete|partial) genome\b|\bmitogenomes?\b|\bmitochondrial genomes?\b', title)
                or (re.search(r'\b(?:mitochondrion|mitochondria)\b', title) and 'complete sequence' in title))


def filter_references(download, output, report, kept_accessions, include_mitogenomes, reference_taxa):
    """Apply the whole-genome title policy, label references with their taxon and audit every decision."""
    kept = excluded = 0
    with Path(output).open('x') as writer, Path(report).open('x', newline='') as audit, Path(kept_accessions).open('x') as ids:
        table = csv.writer(audit, delimiter='\t', lineterminator='\n')
        table.writerow(['accession', 'description', 'length_nt', 'decision', 'reason'])
        for header, sequence in fasta(download):
            identifier = header.split()[0]
            genome = genome_description(header)
            retain = include_mitogenomes or not genome
            table.writerow([identifier, header, len(sequence), 'retained' if retain else 'excluded',
                            'genome_description' if genome else 'no_genome_description'])
            if not retain:
                excluded += 1
                continue
            base = reference_id(identifier)
            accession = next((x for x in [base] + base.split('|') if x in reference_taxa), None)
            if accession is None:
                raise ValueError(f'No taxon metadata for reference {identifier}')
            write_record(writer, labelled_reference_id(identifier, reference_taxa[accession]), sequence)
            ids.write(identifier + '\n')
            kept += 1
    if not kept:
        raise ValueError('No references remain after mitogenome filtering; inspect reference_filter.tsv')
    return kept, excluded


def filter_translations(representatives, translated, output, report, minimum):
    """Keep sequences with at least `minimum` consecutive amino acids free of stops (*) and X in any frame."""
    lengths = {}
    for header, peptide in fasta(translated):
        match = re.fullmatch(r'(.+)_frame=(-?[123])', header.split()[0])
        if not match:
            raise ValueError(f'Unexpected SeqKit translation header: {header}')
        identifier, frame = match.groups()
        frames = lengths.setdefault(identifier, {})
        if frame in frames:
            raise ValueError(f'Duplicate translated frame: {header}')
        longest = max(map(len, re.split(r'[*X]', peptide.upper())), default=0)
        frames[frame] = (len(peptide), longest)
    kept, rows, seen = 0, [], set()
    with Path(output).open('x') as handle:
        for header, sequence in fasta(representatives):
            identifier = header.split()[0]
            if identifier in seen:
                raise ValueError(f'Duplicate representative ID: {identifier}')
            seen.add(identifier)
            frames = lengths.get(identifier, {})
            if set(frames) != {'1', '2', '3', '-1', '-2', '-3'}:
                raise ValueError(f'Missing translation frames for {identifier}')
            passed = any(longest >= minimum for _, longest in frames.values())
            for frame, (total, longest) in frames.items():
                rows.append((identifier, frame, total, longest, minimum, passed))
            if passed:
                write_record(handle, header, sequence)
                kept += 1
    if set(lengths) != seen:
        raise ValueError('Translation IDs do not match representatives')
    tsv(report, ['sequence_id', 'frame', 'translation_length', 'longest_stop_free_length',
                 'minimum_amino_acids', 'sequence_pass'], rows)
    return kept


def vsearch_cluster(runner, source, consensus_sequences, uc, stdout, threads, identity):
    """Build consensuses at the requested identity, restoring seed headers and keeping UC memberships."""
    lengths = [len(sequence) for _, sequence in fasta(source)]
    raw_consensus = Path(str(consensus_sequences) + '.raw')
    runner.run(['vsearch', '--cluster_fast', source, '--id', identity, '--iddef', '2',
                '--strand', 'both', '--qmask', 'none', '--minseqlength', '1',
                '--maxseqlength', max(lengths), '--notrunclabels',
                '--consout', raw_consensus, '--uc', uc, '--threads', threads], stdout)
    source_headers = {header for header, _ in fasta(source)}
    seen = set()
    with Path(consensus_sequences).open('x') as output:
        for header, sequence in fasta(raw_consensus):
            # VSEARCH adds centroid= and ;seqs=N to consensus labels.
            label, separator, size = header.removeprefix('centroid=').rpartition(';seqs=')
            if not header.startswith('centroid=') or not separator or not size.isdigit() or label not in source_headers:
                raise ValueError(f'Unexpected VSEARCH consensus header: {header}')
            if label in seen:
                raise ValueError(f'Duplicate consensus seed label: {label}')
            seen.add(label)
            write_record(output, label, sequence)
    raw_consensus.unlink()
    return len(lengths), len(seen)


def cluster_references(runner, source, output, groups_dir, summary, memberships, threads):
    """Cluster the references of each species at 98% and keep one consensus per cluster."""
    identity = '0.98'
    groups_dir.mkdir()
    groups, handles, seen = {}, OrderedDict(), set()
    try:
        for header, sequence in fasta(source):
            identifier = header.split()[0]
            taxon = identifier.partition(TAXON_SEPARATOR)[2]
            if not taxon:
                raise ValueError(f'Reference lacks a taxonomic label: {header}')
            key = (taxon,)
            if identifier in seen:
                raise ValueError(f'Duplicate sequence ID: {identifier}')
            seen.add(identifier)
            if key not in groups:
                group_id = f'group{len(groups)+1:06d}'
                groups[key] = (group_id, groups_dir/f'{group_id}.input.fasta')
            handle = handles.pop(key, None)
            if handle is None:
                if len(handles) >= 64:
                    handles.popitem(last=False)[1].close()
                handle = groups[key][1].open('a')
            handles[key] = handle
            write_record(handle, header, sequence)
    finally:
        for handle in handles.values():
            handle.close()
    if not groups:
        raise ValueError('No sequences to cluster')
    rows = []
    with output.open('x') as combined, memberships.open('x', newline='') as audit:
        writer = csv.writer(audit, delimiter='\t', lineterminator='\n')
        writer.writerow(['group_id', 'record_type', 'cluster_number', 'size_or_length', 'percent_identity',
                         'strand', 'unused_1', 'unused_2', 'alignment', 'query', 'target'])
        for key, (group_id, path) in groups.items():
            consensus = groups_dir/f'{group_id}.consensus.fasta'
            uc = groups_dir/f'{group_id}.uc'
            before, after = vsearch_cluster(runner, path, consensus, uc, groups_dir/f'{group_id}.stdout.log',
                                            threads, identity)
            for header, sequence in fasta(consensus):
                write_record(combined, header, sequence)
            with uc.open() as handle:
                for line in handle:
                    if line.strip() and not line.startswith('#'):
                        writer.writerow([group_id] + line.rstrip('\n').split('\t'))
            rows.append((group_id, *key, before, after, identity))
    tsv(summary, ['group_id', 'taxon_label', 'input_sequences', 'consensus_sequences', 'identity'], rows)
    runner.note(f'VSEARCH: {len(seen)} sequences in {len(groups)} groups -> {sum(r[-2] for r in rows)} consensus sequences')


def file_digest(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def validate_outputs(paths):
    """Reject missing, empty, malformed or duplicate-ID stage outputs."""
    for path in map(Path, paths):
        if path.is_dir():
            continue
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f'Missing or empty output: {path}')
        if path.suffix in ('.fasta', '.faa'):
            seen = set()
            for header, sequence in fasta(path):
                identifier = header.split()[0]
                if identifier in seen or not sequence:
                    raise ValueError(f'Duplicate ID or empty sequence in {path}: {identifier}')
                seen.add(identifier)
        elif path.name.endswith('.gb.xml'):
            for _ in iter_genbank_records(path):
                pass
        elif path.suffix == '.xml':
            for _, element in ET.iterparse(path, events=['end']):
                element.clear()


def check_alignment(alignment, retained, trimmed=False):
    """Rectangular alignment with the retained IDs; untrimmed rows must equal their input (either strand)."""
    aligned, originals = dict(fasta(alignment)), dict(fasta(retained))
    if len({len(seq) for seq in aligned.values()}) != 1:
        raise ValueError(f'Incomplete alignment: unequal sequence lengths in {alignment}')
    if not set(aligned) <= set(originals) or (not trimmed and set(aligned) != set(originals)):
        raise ValueError(f'Alignment IDs do not match the retained FASTA: {alignment}')
    if any(set(seq.upper()) - set('ACGTRYSWKMBDHVN-') for seq in aligned.values()):
        raise ValueError(f'Invalid nucleotide characters in {alignment}')
    if trimmed:
        return
    complement = str.maketrans('ACGTRYMKBDHVN', 'TGCAYRKMVHDBN')
    for name, seq in aligned.items():
        ungapped, original = seq.replace('-', '').upper(), originals[name].upper()
        if ungapped not in (original, original.translate(complement)[::-1]):
            raise ValueError(f'Alignment sequence differs from retained sequence: {name}')


class Runner:
    def __init__(self, directory):
        self.directory = directory
        self.log = directory/'postprocessing.log'

    def note(self, message):
        print(message, flush=True)
        with self.log.open('a') as handle:
            handle.write(message + '\n')

    def run(self, command, output, input_path=None, cwd=None):
        command = list(map(str, command))
        self.note('$ ' + shlex.join(command))
        source = Path(input_path).open('rb') if input_path else subprocess.DEVNULL
        try:
            with Path(output).open('xb') as target:
                result = subprocess.run(command, stdin=source, stdout=target, stderr=subprocess.PIPE, cwd=cwd)
            with self.log.open('ab') as log:
                log.write(result.stderr)
            if result.returncode:
                # Show the tool's own last error lines, not just the exit status.
                tail = result.stderr.decode(errors='replace').strip().splitlines()[-5:]
                raise ValueError(f'{Path(command[1] if command[0] == sys.executable else command[0]).name} '
                                 f'failed (exit {result.returncode}):\n  ' + '\n  '.join(tail))
        finally:
            if input_path:
                source.close()


class Stages:
    """Run each stage, or on --resume reuse it when outputs, inputs and settings are unchanged."""

    def __init__(self, runner, resume):
        self.runner = runner
        self.path = runner.directory/'checkpoints.json'
        self.records = json.loads(self.path.read_text()) if resume and self.path.exists() else {}

    @staticmethod
    def fingerprint(paths):
        # Size and modification time detect regenerated inputs without re-hashing large files.
        return {str(p): [p.stat().st_size, p.stat().st_mtime_ns] for p in map(Path, paths)}

    def run(self, name, outputs, build, inputs=(), settings=None, extra=()):
        outputs = [Path(p) for p in outputs]
        expected = {'inputs': self.fingerprint(inputs), 'settings': settings}
        record = self.records.get(name)
        if record and all(record.get(key) == value for key, value in expected.items()):
            try:
                validate_outputs(outputs)
                if record['outputs'] == {p.name: file_digest(p) for p in outputs if p.is_file()}:
                    self.runner.note(f'REUSE {name}: outputs, inputs and settings unchanged')
                    return
            except (OSError, ValueError, ET.ParseError):
                pass
        stale = [p for p in dict.fromkeys(outputs + [Path(x) for x in extra]) if p.exists()]
        if stale:
            archive = self.runner.directory/'resume_archive'/f'{time.time_ns()}_{name}'
            archive.mkdir(parents=True)
            for path in stale:
                path.rename(archive/path.name)
            self.runner.note(f'Archived previous {name} outputs in {archive}')
        (self.runner.directory/'postprocessing.complete').unlink(missing_ok=True)
        build()
        validate_outputs(outputs)
        self.records[name] = {**expected, 'outputs': {p.name: file_digest(p) for p in outputs if p.is_file()}}
        temporary = self.path.with_suffix('.tmp')
        temporary.write_text(json.dumps(self.records, indent=2) + '\n')
        temporary.replace(self.path)


def check_environment(args):
    required = ['vsearch', 'seqkit', 'clipkit', 'iqtree'] + (['kc-align', 'kalign', 'mafft'] if args.kcalign else ['mafft'])
    if args.reference_source == 'download':
        required += ['esearch', 'efetch']
    missing = [name for name in required if not shutil.which(name)]
    if missing:
        raise ValueError('Missing programs: ' + ', '.join(missing) + '. Activate the taxonomy environment.')
    from importlib.util import find_spec
    if find_spec('baltic') is None:
        raise ValueError('Baltic is required: conda install -n taxonomy -c bioconda baltic')
    return {name: shutil.which(name) for name in required}


def kcalign_threads(threads):
    """kc-align mixed mode splits work in threes: use 1 or the largest multiple of 3 not above --threads."""
    return 1 if threads < 3 else threads - threads % 3


def align_mafft(runner, retained, alignment, orientation, threads):
    """MAFFT FFT-NS-2 with default gap penalties; restore headers of reverse-complemented rows (_R_) and record their strand."""
    partial = Path(str(alignment) + '.partial')
    partial.unlink(missing_ok=True)
    try:
        runner.run(['mafft', '--retree', '2', '--maxiterate', '0', '--adjustdirectionaccurately',
                    '--thread', threads, retained], partial)
        originals = {header for header, _ in fasta(retained)}
        restored, strands = [], []
        for header, sequence in fasta(partial):
            reverse = header not in originals and header.startswith('_R_') and header[3:] in originals
            original = header[3:] if reverse else header
            if original not in originals:
                raise ValueError(f'Unexpected MAFFT header: {header}')
            restored.append((original, sequence))
            strands.append((original.split()[0], '-' if reverse else '+'))
        with alignment.open('x') as handle:
            for header, sequence in restored:
                write_record(handle, header, sequence)
        tsv(orientation, ['sequence_id', 'alignment_strand'], strands)
    finally:
        partial.unlink(missing_ok=True)


def trim_alignment(runner, alignment, trimmed, log):
    """ClipKIT gappy mode: remove alignment columns that are >=90% gaps. Sequences left with
    no base at all (they only covered removed columns) are dropped, because IQ-TREE refuses them."""
    clipkit_output = Path(str(trimmed) + '.clipkit')
    clipkit_output.unlink(missing_ok=True)
    runner.run(['clipkit', alignment, '-m', 'gappy', '-s', 'nt', '-if', 'fasta', '-of', 'fasta', '-o', clipkit_output], log)
    original, result = dict(fasta(alignment)), dict(fasta(clipkit_output))
    clipkit_output.unlink()
    if set(result) != set(original):
        raise ValueError('ClipKIT changed the sequence set')
    if len({len(seq) for seq in result.values()}) != 1 or not next(iter(result.values())):
        raise ValueError('ClipKIT produced an empty or unequal-length alignment')
    empty = [name for name, seq in result.items() if not set(seq.upper()) & set('ACGT')]
    if empty:
        runner.note(f'Removed {len(empty)} sequences with no bases left after trimming: ' + ', '.join(empty))
    if len(result) - len(empty) < 4:
        raise ValueError('Fewer than four sequences have bases left after trimming')
    with Path(trimmed).open('x') as handle:
        for name, seq in result.items():
            if name not in empty:
                write_record(handle, name, seq)


def resolve_taxon(args, root, stages, runner):
    """BOLD uses the plain name; NCBI resolves the name or taxid in the NCBI taxonomy."""
    taxon = args.taxon.strip()
    if not taxon or any(char in taxon for char in '\n\r"[]'):
        raise ValueError('Taxon must be a numeric taxid or a plain taxonomic name')
    if args.reference_source == 'bold':
        if taxon.isdigit():
            raise ValueError('--bold requires a taxon name; use --download for numeric taxids')
        return taxon, taxon
    taxonomy, search = root/'taxonomy.xml', root/'taxonomy_search.xml'

    def fetch():
        if taxon.isdigit():
            runner.run(['efetch', '-db', 'taxonomy', '-id', taxon, '-format', 'xml'], taxonomy)
        else:
            runner.run(['esearch', '-db', 'taxonomy', '-query', f'"{taxon}"[All Names]'], search)
            runner.run(['efetch', '-format', 'xml'], taxonomy, search)
    stages.run('taxonomy', [taxonomy], fetch, settings=taxon, extra=[search])
    taxa = ET.parse(taxonomy).getroot().findall('./Taxon')
    if len(taxa) != 1:
        raise ValueError(f'Taxon resolved to {len(taxa)} NCBI records; supply an unambiguous taxid')
    taxid, name = taxa[0].findtext('TaxId'), taxa[0].findtext('ScientificName')
    if not taxid or not name:
        raise ValueError('NCBI taxonomy response lacks taxid/name')
    return taxid, name


def run_pipeline(args):
    """Check the inputs, then run or resume every stage."""
    if args.kcalign:
        if args.kcalign_reference is None:
            raise ValueError('--kcalign requires --kcalign-reference with a trusted coding FASTA')
        from kcalign_alignment import validate_reference
        validate_reference(args.kcalign_reference, args.translation_table)
    elif args.kcalign_reference is not None:
        raise ValueError('--kcalign-reference requires --kcalign')
    if not (args.analysis_dir/'analysis.complete').is_file():
        raise ValueError('Provide a completed Read processing analysis directory (--analysis-dir)')
    if args.dataset is not None and args.folders:
        raise ValueError('Use a dataset argument or --folders, not both')
    args.folders = analysis_folders([args.dataset] if args.dataset else args.folders, args.analysis_dir)
    executables = check_environment(args)
    bold_files = database_files(args.bold_db or CODE.parent/'data/bold') if args.reference_source == 'bold' else []
    root = args.output_dir
    if args.resume:
        previous = json.loads((root/'parameters.json').read_text())
        args.date = args.date or previous['date']
        for key in FROZEN:
            old, current = previous.get(key), getattr(args, key)
            equal = (Path(old).resolve() == Path(current).resolve()) if key == 'analysis_dir' else old == current
            if not equal:
                raise ValueError(f'Cannot resume with changed {key}: {old!r} -> {current!r}; use a new output directory')
    else:
        args.date = args.date or date.today().strftime('%Y%m%d')
        root.mkdir(parents=True, exist_ok=False)
    runner = Runner(root)
    stages = Stages(runner, args.resume)
    config = {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()}
    config.update(executables=executables, alignment_algorithm='kc-align' if args.kcalign else 'FFT-NS-2',
                  alignment_trimming='clipkit-gappy-nt', bootstrap_method='none')
    (root/('resume_parameters.json' if args.resume else 'parameters.json')).write_text(json.dumps(config, indent=2) + '\n')
    try:
        build_tree(args, root, runner, stages, bold_files)
        (root/'postprocessing.complete').write_text('complete\n')
        (root/'postprocessing.failed').unlink(missing_ok=True)
    except Exception as error:
        runner.note(f'FAILED: {error}')
        (root/'postprocessing.failed').write_text(str(error) + '\n')
        raise


def build_tree(args, root, runner, stages, bold_files):
    taxid, name = resolve_taxon(args, root, stages, runner)
    slug = re.sub(r'[^A-Za-z0-9_.-]+', '_', name).strip('_.')
    stem = f'{slug}_{args.date}'
    shown = name if taxid == name else f'{name} (NCBI taxid {taxid})'
    runner.note(f'Resolved taxon: {shown}; folders: {", ".join(args.folders)}')
    tables = [args.analysis_dir/folder/f'{folder}_clusters.tsv' for folder in args.folders]

    # 1. Local clusters (consensus sequences and single reads from Module 1) confidently
    #    assigned to the taxon, then the main identification of each specimen.
    local, extraction_log = root/f'clusters_{stem}.fasta', root/'extraction.log'
    stages.run('extraction', [local, extraction_log], lambda: runner.run(
        [sys.executable, CODE/'export_clusters.py', '--taxon', taxid, '--taxon', name,
         '--folders', *args.folders, '--analysis-dir', args.analysis_dir, '--output', local], extraction_log),
        inputs=tables, settings=sorted({taxid, name}))
    runner.note(extraction_log.read_text().rstrip())
    summary, local_reps = root/'extraction_by_specimen_taxon.tsv', root/f'clusters_{stem}.main_taxon.fasta'

    def main_taxon():
        runner.note(f'Main identification per specimen: kept {select_main_taxon(local, summary, local_reps)} clusters')
    stages.run('main_taxon_selection', [summary, local_reps], main_taxon, inputs=[local])

    # 2. External references, labelled with their taxon, length-filtered and clustered per taxon.
    prefix_name = 'BOLD' if args.reference_source == 'bold' else 'COI'
    accessions = root/'reference_accessions.txt'
    if args.reference_source == 'bold':
        download, metadata = root/'references_bold.fasta', root/'bold_reference_metadata.tsv'
        (root/'reference_query.txt').write_text(f'BOLD header contains (case-insensitive): {name}\n')
        stages.run('bold_extraction', [download, accessions, metadata],
                   lambda: runner.note(f'Extracted {extract_bold(bold_files, name, download, metadata, accessions)} '
                                       f'BOLD records matching {name}'), inputs=bold_files, settings=name)
    else:
        query = reference_query(taxid, args.include_mitogenomes)
        (root/'reference_query.txt').write_text(query + '\n')
        download, metadata, search = root/'references_downloaded.fasta', root/'reference_records.gb.xml', root/'coi_search.xml'

        def fetch():
            runner.run(['esearch', '-db', 'nuccore', '-query', query], search)
            runner.run(['efetch', '-format', 'acc'], accessions, search)
            if not accessions.read_text().strip():
                raise ValueError('NCBI returned no COI accessions')
            runner.run(['efetch', '-format', 'fasta'], download, search)
            runner.run(['efetch', '-db', 'nuccore', '-format', 'gb', '-mode', 'xml'], metadata, accessions)
        stages.run('download', [download, accessions, metadata], fetch, settings=query, extra=[search])
    references = root/f'{prefix_name}_{stem}.fasta'
    policy_outputs = [references, root/'reference_filter.tsv', root/'reference_accessions_retained.txt']

    def policy():
        kept, excluded = filter_references(download, *policy_outputs, args.include_mitogenomes, reference_names(metadata))
        runner.note(f'Reference policy: include_mitogenomes={args.include_mitogenomes}; retained {kept}, excluded {excluded}')
    stages.run('reference_policy', policy_outputs, policy, inputs=[download, metadata], settings=args.include_mitogenomes)
    length_filtered = root/f'{prefix_name}_{stem}.length500-2000.fasta'
    stages.run('reference_lengths', [length_filtered], lambda: runner.run(
        ['seqkit', 'seq', '--min-len', '500', '--max-len', '2000', references], length_filtered), inputs=[references])
    reference_count = sum(1 for _ in fasta(length_filtered))
    if not reference_count:
        raise ValueError('No reference sequences remain in the inclusive 500–2000 nt range')
    runner.note(f'References retained after 500–2000 nt filtering: {reference_count}')
    reference_reps = root/f'{prefix_name}_{stem}.vsearch98.fasta'
    reference_outputs = [reference_reps, root/'vsearch_reference_memberships.tsv', root/'vsearch_reference_groups.tsv']
    stages.run('vsearch_references', reference_outputs, lambda: cluster_references(
        runner, length_filtered, reference_reps, root/'vsearch_reference_groups', reference_outputs[2],
        reference_outputs[1], args.threads), inputs=[length_filtered], extra=[root/'vsearch_reference_groups'])

    # 3. Combine local clusters and reference consensus sequences; keep those with a long open reading frame.
    representatives = root/f'combined_{stem}.fasta'

    def combine():
        with representatives.open('x') as handle:
            for source in (local_reps, reference_reps):
                for header, sequence in fasta(source):
                    write_record(handle, header, sequence)
    stages.run('combined_representatives', [representatives], combine, inputs=[local_reps, reference_reps])
    stage_label = 'selected' if args.skip_translation else 'nonumt'
    prefix = root/f'combined_{slug}_{stage_label}_{args.date}'
    retained = Path(str(prefix) + '.fasta')
    if args.skip_translation:
        stages.run('selection', [retained], lambda: shutil.copyfile(representatives, retained), inputs=[representatives])
        runner.note('Translation and amino-acid selection skipped; retained all consensuses')
    else:
        translated = root/f'combined_{stem}.six_frames.faa'
        stages.run('translation', [translated], lambda: runner.run(
            ['seqkit', 'translate', '-f', '6', '-F', '-T', args.translation_table, representatives], translated),
            inputs=[representatives], settings=args.translation_table)
        report = root/'translation_lengths.tsv'
        stages.run('retained', [retained, report], lambda: filter_translations(
            representatives, translated, retained, report, args.min_aa),
            inputs=[representatives, translated], settings=[args.min_aa, 'stop-free'])
    count = sum(1 for _ in fasta(retained))
    runner.note(f'Retained {count} nucleotide consensus sequences for alignment')
    if count < 4:
        raise ValueError('Fewer than four representatives passed; FASTA saved, phylogeny not run')

    # 4. Align, trim gappy columns, infer the maximum-likelihood tree and draw it.
    alignment = Path(str(prefix) + '.aln.fasta')
    if args.kcalign:
        from kcalign_alignment import align_kcalign
        artifacts = [Path(str(prefix) + s) for s in ('.pre_align.aa.fasta', '.kcalign_selection.tsv', '.kcalign_ids.tsv')]
        stages.run('alignment', [alignment, *artifacts], lambda: align_kcalign(
            runner, retained, alignment, args.kcalign_reference, args.translation_table, kcalign_threads(args.threads)),
            inputs=[retained, args.kcalign_reference], settings=['kc-align', args.translation_table],
            extra=[Path(str(prefix) + '.kcalign.stdout.log'), root/'alignment_orientation.tsv'])
        check_alignment(alignment, retained, trimmed=True)
    else:
        orientation = root/'alignment_orientation.tsv'
        stages.run('alignment', [alignment, orientation], lambda: align_mafft(
            runner, retained, alignment, orientation, args.threads), inputs=[retained], settings=['FFT-NS-2', 'default gap penalties'],
            extra=[Path(str(prefix) + s) for s in ('.pre_align.aa.fasta', '.kcalign_selection.tsv', '.kcalign_ids.tsv',
                                                   '.kcalign.stdout.log')])
        check_alignment(alignment, retained)
    trimmed = Path(str(prefix) + '.clipkit.aln.fasta')
    stages.run('clipkit', [trimmed], lambda: trim_alignment(runner, alignment, trimmed, root/'clipkit.stdout.log'),
               inputs=[alignment], settings=['gappy', 'nt', 'drop empty rows'], extra=[root/'clipkit.stdout.log'])
    treefile = Path(str(prefix) + '.treefile')
    iqtree_files = [Path(str(prefix) + s) for s in ('.iqtree', '.log', '.ckp.gz', '.bionj', '.mldist', '.model.gz', '.uniqueseq.phy')]
    stages.run('phylogeny', [treefile], lambda: runner.run(
        ['iqtree', '-s', trimmed, '-m', 'MFP', '-T', args.threads, '--seed', '666', '--prefix', prefix],
        root/'iqtree.stdout.log'), inputs=[trimmed], settings=['MFP', 666, 'no bootstrap'],
        extra=[root/'iqtree.stdout.log', *iqtree_files])
    figure = [Path(str(prefix) + '.figure' + s) for s in ('.pdf', '.svg', '.png', '.midpoint.nwk', '.tips.tsv', '.plot.json')]
    stages.run('figure', figure, lambda: runner.run(
        [sys.executable, CODE/'plot_phylogeny.py', '--tree', treefile,
         '--reference-metadata', metadata, '--output-prefix', str(prefix) + '.figure'], root/'tree_figure.stdout.log'),
        inputs=[treefile, metadata], extra=[root/'tree_figure.stdout.log'])
    runner.note(f'Completed: {treefile}')


def main():
    """Parse source, filtering and alignment choices before starting or resuming a run."""
    parser = argparse.ArgumentParser(description=__doc__)
    reference_group = parser.add_mutually_exclusive_group(required=True)
    reference_group.add_argument('--download', dest='reference_source', action='store_const', const='download',
                                 help='Download NCBI COI references using Entrez Direct')
    reference_group.add_argument('--bold', dest='reference_source', action='store_const', const='bold',
                                 help='Select local BOLD FASTA records whose header contains the taxon name; no network')
    parser.add_argument('--bold-db', type=Path, help='BOLD gzipped FASTA file or directory (default: data/bold)')
    parser.add_argument('--include-mitogenomes', action='store_true', help='Keep whole-mitogenome records (excluded by default)')
    parser.add_argument('--taxon', required=True, help='Taxon name (or NCBI taxid with --download)')
    parser.add_argument('dataset', nargs='?', help='Dataset parent, e.g. sequencing2: use all its analysed folders')
    parser.add_argument('--folders', nargs='+', help='Sample folders; default: all folders in the completed analysis')
    parser.add_argument('--analysis-dir', type=Path, required=True, help='Read processing analysis directory')
    parser.add_argument('--output-dir', type=Path, required=True, help='New output directory, or an existing run with --resume')
    parser.add_argument('--resume', action='store_true',
                        help='Continue an existing run: reuse unchanged stages, redo stages whose inputs or settings changed')
    parser.add_argument('--date', default=None, help='YYYYMMDD used in file names; default today, or the saved date with --resume')
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--kcalign-reference', type=Path, help='Trusted single-sequence, in-frame coding FASTA required with --kcalign')
    parser.add_argument('--skip-translation', action='store_true',
                        help='Skip the six-frame translation screen; keep all consensuses (e.g. non-coding loci)')
    aligners = parser.add_mutually_exclusive_group()
    aligners.add_argument('--kcalign', action='store_true', help='Use kc-align codon alignment instead of MAFFT')
    aligners.add_argument('--mafft', dest='kcalign', action='store_false', help='Use MAFFT FFT-NS-2 (default)')
    parser.add_argument('--translation-table', type=int, default=5, help='Genetic code (5: invertebrate mitochondrial)')
    parser.add_argument('--min-aa', type=int, default=200,
                        help='Minimum number of consecutive amino acids without a stop codon or X, in any frame')
    args = parser.parse_args()
    if args.bold_db is not None and args.reference_source != 'bold':
        parser.error('--bold-db requires --bold')
    if min(args.threads, args.min_aa, args.translation_table) < 1:
        parser.error('Numeric settings must be positive')
    if args.date is not None and not re.fullmatch(r'\d{8}', args.date):
        parser.error('--date must have YYYYMMDD format')
    if args.resume and not (args.output_dir/'parameters.json').is_file():
        parser.error('--resume requires an existing run directory with parameters.json')
    try:
        run_pipeline(args)
    except (ValueError, OSError, subprocess.CalledProcessError, ET.ParseError) as error:
        parser.exit(1, f'Post-processing failed: {error}\n')


if __name__ == '__main__':
    main()
