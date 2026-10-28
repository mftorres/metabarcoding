#!/usr/bin/env python3
"""Export a SINTAX taxon, retrieve COI references, cluster, filter and infer a tree."""
import argparse
from collections import Counter, OrderedDict
import csv
import hashlib
from contextlib import contextmanager
import time
import tempfile
from datetime import date
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
from urllib.parse import unquote
import xml.etree.ElementTree as ET
from genbank_xml import iter_genbank_records
from sequence_headers import HEADER_SCHEMA
from bold_references import database_files, extract_bold
from export_sintax_fasta import TAXON_SELECTION_POLICY


def fasta(path):
    """Stream multiline FASTA while preserving complete tracking headers."""
    header, parts = None, []
    with Path(path).open() as handle:
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


def write_record(handle, header, sequence):
    """Write a FASTA record with consistent 80-base wrapping."""
    handle.write(f'>{header}\n')
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start:start + 80] + '\n')


def tsv(path, columns, rows):
    """Write a tab-delimited audit table with its column headers."""
    with Path(path).open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(columns)
        writer.writerows(rows)


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
                or (re.search(r'\b(?:mitochondrion|mitochondria)\b', title)
                    and 'complete sequence' in title))


def filter_references(download, output, report, kept_accessions, include_mitogenomes=False):
    """Apply the reference title policy and record every keep/exclude decision."""
    kept = excluded = 0
    with Path(output).open('x') as writer, Path(report).open('x', newline='') as audit, Path(kept_accessions).open('x') as ids:
        table = csv.writer(audit, delimiter='\t')
        table.writerow(['accession', 'description', 'length_nt', 'decision', 'reason'])
        for header, sequence in fasta(download):
            genome = genome_description(header)
            retain = include_mitogenomes or not genome
            table.writerow([header.split()[0], header, len(sequence), 'retained' if retain else 'excluded',
                            'genome_description' if genome else 'no_genome_description'])
            if retain:
                write_record(writer, header, sequence)
                ids.write(header.split()[0] + '\n')
                kept += 1
            else:
                excluded += 1
    if not kept:
        raise ValueError('No references remain after mitogenome filtering; inspect reference_filter.tsv')
    return kept, excluded


def filter_translations(representatives, translated, output, report, minimum, mode):
    """Measure all six frames and retain original nucleotides meeting the peptide rule."""
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
            passed = any(values[0 if mode == 'total' else 1] >= minimum for values in frames.values())
            for frame, (total, longest) in frames.items():
                rows.append((identifier, frame, total, longest, mode, minimum, passed))
            if passed:
                write_record(handle, header, sequence)
                kept += 1
    if set(lengths) != seen:
        raise ValueError('Translation IDs do not match representatives')
    tsv(report, ['sequence_id', 'frame', 'translation_length', 'longest_stop_free_length',
                 'length_mode', 'minimum_amino_acids', 'sequence_pass'], rows)
    return kept


CLUSTERING_POLICY = 'vsearch98-consensus-local-barcode-forward-reverse-taxon-reference-global-v3'


def vsearch_cluster(runner, source, consensus_sequences, uc, stdout, threads):
    """Build 98% consensuses and restore seed headers while retaining UC memberships."""
    lengths = [len(sequence) for _, sequence in fasta(source)]
    if not lengths:
        raise ValueError(f'No sequences to cluster in {source}')
    raw_consensus = Path(str(consensus_sequences) + '.raw')
    runner.run(['vsearch', '--cluster_fast', source, '--id', '0.98', '--iddef', '2',
                '--strand', 'both', '--qmask', 'none', '--minseqlength', '1',
                '--maxseqlength', max(lengths), '--notrunclabels',
                '--consout', raw_consensus, '--uc', uc, '--threads', threads], stdout)
    validate_outputs([raw_consensus, uc])
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
    validate_outputs([consensus_sequences, uc])
    return len(lengths), sum(1 for _ in fasta(consensus_sequences))


def cluster_local_groups(runner, source, output, groups_dir, summary, memberships, threads):
    """Keep barcode/forward/reverse/taxon groups separate throughout local clustering."""
    groups_dir.mkdir()
    groups, handles, seen = {}, OrderedDict(), set()
    try:
        for header, sequence in fasta(source):
            identifier = header.split()[0]
            parts = identifier.split('-', 5)
            if len(parts) != 6 or parts[4] not in ('contamination', 'nocontamination'):
                raise ValueError(f'Expected a primer-aware local FASTA header: {header}')
            if identifier in seen:
                raise ValueError(f'Duplicate local sequence ID: {identifier}')
            seen.add(identifier)
            key = tuple(parts[:4])
            if key not in groups:
                group_id = f'group{len(groups)+1:06d}'
                groups[key] = (group_id, groups_dir / f'{group_id}.input.fasta')
            group_id, path = groups[key]
            handle = handles.pop(key, None)
            if handle is None:
                if len(handles) >= 64:
                    _, old = handles.popitem(last=False); old.close()
                handle = path.open('a')
            handles[key] = handle
            write_record(handle, header, sequence)
    finally:
        for handle in handles.values():
            handle.close()
    if not groups:
        raise ValueError('No local reads to cluster')
    rows = []
    with output.open('x') as combined, memberships.open('x', newline='') as audit:
        writer = csv.writer(audit, delimiter='\t')
        writer.writerow(['group_id', 'record_type', 'cluster_number', 'size_or_length',
                         'percent_identity', 'strand', 'unused_1', 'unused_2', 'alignment', 'query', 'target'])
        for key, (group_id, path) in groups.items():
            consensus_sequences = groups_dir / f'{group_id}.consensus.fasta'
            uc = groups_dir / f'{group_id}.uc'
            before, after = vsearch_cluster(runner, path, consensus_sequences, uc,
                                           groups_dir / f'{group_id}.stdout.log', threads)
            for header, sequence in fasta(consensus_sequences):
                if tuple(header.split('-', 4)[:4]) != key:
                    raise ValueError(f'VSEARCH changed grouping fields: {header}')
                write_record(combined, header, sequence)
            with uc.open() as handle:
                for line in handle:
                    if line.strip() and not line.startswith('#'):
                        writer.writerow([group_id] + line.rstrip('\n').split('\t'))
            rows.append((group_id, *(unquote(value) for value in key), before, after, '0.98'))
    tsv(summary, ['group_id', 'source_folder', 'forward_primer', 'reverse_primer', 'taxon_label',
                  'input_sequences', 'consensus_sequences', 'identity'], rows)
    runner.note(f'Local VSEARCH: {len(seen)} reads in {len(groups)} barcode/primer/taxon groups -> '
                f'{sum(row[6] for row in rows)} consensus sequences')


def file_digest(path):
    """Hash output contents for resume validation without loading the whole file."""
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def validate_outputs(paths):
    """Reject missing, empty, malformed or duplicate-ID stage outputs."""
    for path in paths:
        path = Path(path)
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f'Missing or empty output: {path}')
        if path.suffix in ('.fasta', '.faa'):
            seen = set()
            for header, sequence in fasta(path):
                identifier = header.split()[0]
                if identifier in seen or not sequence:
                    raise ValueError(f'Duplicate ID or empty sequence in {path}: {identifier}')
                seen.add(identifier)
            if not seen:
                raise ValueError(f'No FASTA records in {path}')
        elif path.name.endswith('.gb.xml'):
            for _ in iter_genbank_records(path):
                pass
        elif path.suffix == '.xml':
            # Streaming parse also detects truncated downloads without retaining sequences.
            for _, element in ET.iterparse(path, events=['end']):
                element.clear()


def validate_retained_checkpoint(representatives, translated, retained, report, minimum, mode):
    # Older runs have no hashes: independently verify the retained file and report
    # against all six saved translations before adopting them as a checkpoint.
    """Recompute translation selection to validate uncheckpointed saved results."""
    with tempfile.TemporaryDirectory(prefix='sintax_resume_') as directory:
        expected = Path(directory) / 'retained.fasta'
        expected_report = Path(directory) / 'translation_lengths.tsv'
        filter_translations(representatives, translated, expected, expected_report, minimum, mode)
        if file_digest(expected) != file_digest(retained) or file_digest(expected_report) != file_digest(report):
            raise ValueError('Retained sequences/report do not match saved translations')


def valid_alignment(alignment, retained=None, allow_trimmed=False):
    """Check rectangular FASTA and IDs; MAFFT must also preserve nucleotide content."""
    validate_outputs([alignment])
    aligned = dict(fasta(alignment))
    if len({len(seq) for seq in aligned.values()}) != 1:
        raise ValueError('Incomplete alignment: unequal sequence lengths')
    if retained is not None and retained.exists():
        originals = dict(fasta(retained))
        ids_match = set(aligned).issubset(originals) if allow_trimmed else set(aligned) == set(originals)
        if not ids_match:
            raise ValueError('Incomplete alignment: sequence IDs do not match retained FASTA')
        if allow_trimmed:
            if any(set(seq.upper()) - set('ACGTRYSWKMBDHVN-') for seq in aligned.values()):
                raise ValueError('Invalid nucleotide characters in trimmed alignment')
            return
        complement = str.maketrans('ACGTRYMKBDHVNacgtrymkbdhvn', 'TGCAYRKMVHDBNtgcayrkmvhdbn')
        for name, seq in aligned.items():
            ungapped = seq.replace('-', '').upper()
            original = originals[name].upper()
            if ungapped not in (original, original.translate(complement)[::-1]):
                raise ValueError(f'Alignment sequence differs from retained sequence: {name}')


def latest_sequences(root, args, records):
    """Read-only discovery of the latest usable sequence stage, without network."""
    stage_label = 'selected' if getattr(args, 'skip_translation', False) else 'nonumt'
    candidates = list(root.glob(f'combined_*_{stage_label}_{args.date}.fasta'))
    # An alignment can also stand alone after earlier intermediate files were removed.
    candidates += [Path(str(p).replace('.aln.fasta', '.fasta'))
                   for p in root.glob(f'combined_*_{stage_label}_{args.date}.aln.fasta') if
                   Path(str(p).replace('.aln.fasta', '.fasta')) not in candidates]
    if len(candidates) > 1:
        raise ValueError('Multiple retained sequence sets found; use a separate run directory per taxon')
    if not candidates:
        return None
    retained = candidates[0]
    slug = retained.name[len('combined_'):-len(f'_{stage_label}_{args.date}.fasta')]
    alignment = retained.with_suffix('.aln.fasta')
    orientation = root / 'alignment_orientation.tsv'
    marker = root / 'alignment.complete.json'
    try:
        invalidated = root / 'invalidated_stages.json'
        if invalidated.exists() and 'alignment' in json.loads(invalidated.read_text()):
            raise ValueError('Alignment was invalidated by an explicit restart')
        valid_alignment(alignment, retained, allow_trimmed=getattr(args, 'kcalign', False))
        if marker.exists():
            saved = json.loads(marker.read_text())
            if saved.get('algorithm') != ('kc-align' if getattr(args, 'kcalign', False) else 'FFT-NS-2'):
                raise ValueError('Alignment algorithm mismatch')
            if retained.exists() and saved.get('retained_sha256') != file_digest(retained):
                raise ValueError('Alignment input hash mismatch')
            if saved['alignment_sha256'] != file_digest(alignment):
                raise ValueError('Alignment completion hash mismatch')
        elif getattr(args, 'kcalign', False):
            raise ValueError('kc-align requires an alignment completion marker')
        elif 'alignment' in records:
            if records['alignment'] != {p.name: file_digest(p) for p in (alignment, orientation)}:
                raise ValueError('Alignment checkpoint hash mismatch')
        elif not orientation.is_file():
            raise ValueError('No proof that alignment finished successfully')
        return slug, retained, True
    except (OSError, ValueError, KeyError):
        pass
    try:
        validate_outputs([retained])
        if getattr(args, 'skip_translation', False):
            record = records.get('selection')
            if record != {retained.name: file_digest(retained)}:
                raise ValueError('Unfiltered selection checkpoint mismatch')
            return slug, retained, False
        report = root / 'translation_lengths.tsv'
        record = records.get('retained')
        if record and record != {p.name: file_digest(p) for p in (retained, report)}:
            raise ValueError('Retained sequence checkpoint mismatch')
        # Validate legacy retained files against their report/translations when available.
        reps = root / f'combined_{slug}_{args.date}.fasta'
        translations = root / f'combined_{slug}_{args.date}.six_frames.faa'
        if not record and all(p.is_file() for p in (report, reps, translations)):
            validate_retained_checkpoint(reps, translations, retained, report, args.min_aa, args.aa_length_mode)
        return slug, retained, False
    except (OSError, ValueError):
        return None


class Checkpoints:
    def __init__(self, runner, resume):
        self.runner, self.resume = runner, resume
        self.path = runner.directory / 'checkpoints.json'
        self.records = json.loads(self.path.read_text()) if self.path.exists() else {}
        self.legacy = resume and not self.path.exists()
        self.dirty = False
        self.invalidated_path = runner.directory / "invalidated_stages.json"
        self.invalidated = set(json.loads(self.invalidated_path.read_text())) if self.invalidated_path.exists() else set()

    @contextmanager
    def stage(self, name, paths, extra=(), legacy=True, legacy_validator=None):
        paths = list(map(Path, paths))
        record = self.records.get(name)
        reusable = False
        independent = name in ('download', 'bold_extraction', 'reference_metadata')
        if self.resume and (not self.dirty or independent) and (name not in self.invalidated or independent) and (record or legacy):
            try:
                validate_outputs(paths)
                if not record and legacy_validator is not None:
                    legacy_validator()
                fingerprints = {str(p.name): file_digest(p) for p in paths}
                reusable = record == fingerprints if record else True
            except (OSError, ValueError, ET.ParseError):
                reusable = False
        if reusable:
            self.runner.note(f'REUSE {name}: validated existing outputs')
        else:
            # Invalidate downstream stages after rebuilding an upstream stage.
            self.dirty = True
            # Persist invalidation before doing work, including interrupted retries.
            if name in self.records:
                keys = list(self.records)
                self.invalidated.update(keys[keys.index(name):])
                self.invalidated_path.write_text(json.dumps(sorted(self.invalidated)))
                self.records = {key: self.records[key] for key in keys[:keys.index(name)]}
            self.path.write_text(json.dumps(self.records, indent=2) + '\n')
            (self.runner.directory / 'postprocessing.complete').unlink(missing_ok=True)
            existing = [p for p in dict.fromkeys(paths + list(map(Path, extra))) if p.exists()]
            if existing:
                archive = self.runner.directory / 'resume_archive' / f'{time.time_ns()}_{name}'
                archive.mkdir(parents=True)
                for path in existing:
                    path.rename(archive / path.name)
                self.runner.note(f'Archived incomplete/stale {name} outputs in {archive}')
        yield not reusable
        validate_outputs(paths)
        self.records[name] = {str(p.name): file_digest(p) for p in paths}
        temporary = self.path.with_suffix('.tmp')
        temporary.write_text(json.dumps(self.records, indent=2) + '\n')
        temporary.replace(self.path)
        self.invalidated.discard(name)
        self.invalidated_path.write_text(json.dumps(sorted(self.invalidated)))


class Runner:
    def __init__(self, directory):
        self.directory = directory
        self.log = directory / 'postprocessing.log'

    def note(self, message):
        print(message, flush=True)
        with self.log.open('a') as handle:
            handle.write(message + '\n')

    def run(self, command, output, input_path=None, cwd=None):
        command = list(map(str, command))
        self.note('$ ' + shlex.join(command))
        source = Path(input_path).open('rb') if input_path else subprocess.DEVNULL
        try:
            with Path(output).open('xb') as target, self.log.open('ab') as log:
                subprocess.run(command, stdin=source, stdout=target, stderr=log, check=True, cwd=cwd)
        finally:
            if input_path:
                source.close()


def run_pipeline(args):
    # Freeze biological settings on resume so cached data cannot bypass new choices.
    """Orchestrate source selection, consensus construction, alignment, tree and figure."""
    args.restart_alignment = getattr(args, 'restart_alignment', False)
    if args.restart_alignment and not getattr(args, 'resume', False):
        raise ValueError('--restart-alignment requires --resume')
    args.skip_translation = getattr(args, 'skip_translation', False)
    args.kcalign = getattr(args, 'kcalign', False)
    args.kcalign_reference = getattr(args, 'kcalign_reference', None)
    reference_snapshot = None
    if args.kcalign:
        if args.kcalign_reference is None:
            raise ValueError('--kcalign requires --kcalign-reference with a trusted coding FASTA')
        if args.threads != 1 and args.threads % 3:
            raise ValueError('kc-align mixed mode requires --threads 1 or a multiple of 3 (e.g. 6)')
        from kcalign_alignment import validate_reference
        validate_reference(args.kcalign_reference, args.translation_table)
        reference_snapshot = file_digest(args.kcalign_reference)
    elif args.kcalign_reference is not None:
        raise ValueError('--kcalign-reference requires --kcalign')
    stage_label = 'selected' if args.skip_translation else 'nonumt'
    args.reference_source = getattr(args, 'reference_source', 'download')
    args.bold_db = getattr(args, 'bold_db', None)
    if args.reference_source == 'bold' and args.taxon.strip().isdigit():
        raise ValueError('--bold requires a taxon name for header matching; use --download for numeric taxids')
    bold_files = database_files(args.bold_db or Path(__file__).resolve().parents[1] / 'data/bold') if args.reference_source == 'bold' else []
    bold_snapshot = [{'path': str(p.resolve()), 'bytes': p.stat().st_size, 'mtime_ns': p.stat().st_mtime_ns} for p in bold_files]
    required = ['vsearch', 'seqkit', 'clipkit', 'kc-align' if args.kcalign else 'mafft', 'iqtree']
    if args.reference_source == 'download':
        required += ['esearch', 'efetch']
    if args.kcalign:
        required += ['kalign', 'mafft']
    missing = [name for name in required if not shutil.which(name)]
    if missing:
        raise ValueError('Missing executables on PATH: ' + ', '.join(missing) +
                         '. Activate taxonomy; install missing packages (see docs/post_processing.md).')
    from importlib.util import find_spec
    from importlib.metadata import version
    if find_spec('baltic') is None or int(version('baltic').split('.')[0]) < 1:
        raise ValueError('Baltic >=1.0.0 is required; conda install -n taxonomy -c conda-forge -c bioconda baltic')
    if not (args.analysis_dir / 'analysis.complete').is_file():
        raise ValueError('Provide a completed SINTAX analysis directory')
    folders = [Path(folder).name for folder in args.folders]
    if len(folders) != len(set(folders)):
        raise ValueError('Duplicate folder names')
    resume = getattr(args, 'resume', False)
    if resume:
        if not (args.output_dir / 'parameters.json').is_file():
            raise ValueError('--resume requires an existing run with parameters.json')
        previous = json.loads((args.output_dir / 'parameters.json').read_text())
        if previous.get('reference_source', 'download') != args.reference_source:
            raise ValueError('Cannot change --download/--bold during resume; use a new output directory')
        if args.reference_source == 'bold' and previous.get('bold_database_files') != bold_snapshot:
            raise ValueError('BOLD database files changed since this run; use a new output directory')
        if not args.restart_alignment and previous.get('kcalign_reference_sha256') != reference_snapshot:
            raise ValueError('Cannot resume with changed trusted alignment reference; use a new output directory')
        if previous.get('clustering_policy') != CLUSTERING_POLICY:
            raise ValueError('This run uses an older clustering policy. Use a new output directory for grouped VSEARCH 98% consensus sequences.')
        if previous.get('taxon_selection_policy') != TAXON_SELECTION_POLICY:
            raise ValueError('This run used an older taxon-selection policy. Use a new output directory to require queried-label confidence >=0.80.')
        if previous.get('header_schema') != HEADER_SCHEMA:
            raise ValueError('This run uses old FASTA headers without primer fields. Use a new output directory to regenerate exports, clusters, alignment and tree with the current header format.')
        if args.date is None:
            args.date = previous['date']
        for key in ('taxon', 'folders', 'analysis_dir', 'date', 'include_mitogenomes',
                    'translation_table', 'min_aa', 'aa_length_mode', 'bootstrap_method', 'skip_translation', 'kcalign'):
            if key == 'kcalign' and args.restart_alignment:
                continue
            current = getattr(args, key)
            old = previous.get(key)
            if key == 'analysis_dir':
                equal = Path(current).resolve() == Path(old).resolve()
            elif key == 'folders':
                equal = sorted(Path(x).name for x in current) == sorted(Path(x).name for x in old)
            else:
                equal = current == old
            if not equal:
                raise ValueError(f'Cannot resume with changed {key}: {old!r} -> {current!r}; use a new output directory')
    else:
        args.date = args.date or date.today().strftime('%Y%m%d')
        args.output_dir.mkdir(parents=True, exist_ok=False)
    runner = Runner(args.output_dir)
    root = args.output_dir
    checkpoints = Checkpoints(runner, resume)
    try:
        config = {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()}
        config['bold_database_files'] = bold_snapshot
        config['executables'] = {name: shutil.which(name) for name in required}
        config['kcalign_reference_sha256'] = reference_snapshot
        config['alignment_trimming'] = 'clipkit-smart-gap-nt-v1'
        config['alignment_algorithm'] = 'kc-align' if args.kcalign else 'FFT-NS-2'
        config['clustering_policy'] = CLUSTERING_POLICY
        config['header_schema'] = HEADER_SCHEMA
        config['taxon_selection_policy'] = TAXON_SELECTION_POLICY
        if not resume:
            (root / 'parameters.json').write_text(json.dumps(config, indent=2) + '\n')
        else:
            (root / 'resume_parameters.json').write_text(json.dumps(config, indent=2) + '\n')
        restart = latest_sequences(root, args, checkpoints.records) if resume else None
        if args.restart_alignment:
            if restart is None or not restart[1].is_file():
                raise ValueError('Restarting alignment requires valid retained/selected consensus sequences')
            # Reuse validated upstream sequences, but invalidate all alignment-dependent results.
            archive = root / 'resume_archive' / f'{time.time_ns()}_alignment_settings'
            archive.mkdir(parents=True)
            shutil.copyfile(root / 'parameters.json', archive / 'parameters.json')
            (root / 'parameters.json').write_text(json.dumps(config, indent=2) + '\n')
            checkpoints.dirty = True
            checkpoints.invalidated.update(('alignment', 'clipkit', 'phylogeny', 'figure'))
            checkpoints.invalidated_path.write_text(json.dumps(sorted(checkpoints.invalidated)))
            (root / 'postprocessing.complete').unlink(missing_ok=True)
            restart = (restart[0], restart[1], False)
        if restart is not None:
            slug, retained, alignment_ready = restart
            stem = f'{slug}_{args.date}'
            runner.note('RESUME directly from completed alignment' if alignment_ready else
                        'RESUME directly from retained/selected consensus sequences')
        else:
            alignment_ready = False
            taxon = args.taxon.strip()
            if not taxon or any(char in taxon for char in '\n\r"[]'):
                raise ValueError('Taxon must be a numeric taxid or a plain taxonomic name')
            if args.reference_source == 'bold':
                name, taxid = taxon, taxon
            else:
                taxonomy = root / 'taxonomy.xml'
                with checkpoints.stage('taxonomy', [taxonomy], extra=[root / 'taxonomy_search.xml'], legacy=True) as execute:
                    if execute:
                        if taxon.isdigit():
                            runner.run(['efetch', '-db', 'taxonomy', '-id', taxon, '-format', 'xml'], taxonomy)
                        else:
                            search = root / 'taxonomy_search.xml'
                            runner.run(['esearch', '-db', 'taxonomy', '-query', f'"{taxon}"[All Names]'], search)
                            runner.run(['efetch', '-format', 'xml'], taxonomy, search)
                taxa = ET.parse(taxonomy).getroot().findall('./Taxon')
                if len(taxa) != 1:
                    raise ValueError(f'Taxon resolved to {len(taxa)} NCBI records; supply an unambiguous taxid')
                taxid, name = taxa[0].findtext('TaxId'), taxa[0].findtext('ScientificName')
                if not taxid or not name:
                    raise ValueError('NCBI taxonomy response lacks taxid/name')
            slug = re.sub(r'[^A-Za-z0-9_.-]+', '_', name).strip('_.')
            stem = f'{slug}_{args.date}'
            runner.note(f'Resolved taxon: {name} (taxid {taxid}); folders: {", ".join(folders)}')
            # Export exact cleaned query IDs using trusted queried-rank SINTAX evidence.
            local = root / f'cleaned_{stem}.fasta'
            extraction_log = root / 'extraction.log'
            with checkpoints.stage('extraction', [local, extraction_log], extra=(), legacy=True) as execute:
                if execute:
                    runner.run([sys.executable, Path(__file__).with_name('export_sintax_taxa.py'),
                                '--taxon', taxid, '--taxon', name, '--folders', *args.folders,
                                '--analysis-dir', args.analysis_dir, '--output', local], extraction_log)
            runner.note(extraction_log.read_text().rstrip())
            counts = Counter(unquote(header.split('-', 1)[0]) for header, _ in fasta(local))
            tsv(root / 'extraction_counts.tsv', ['source_folder', 'reads_extracted'],
                [(folder, counts[folder]) for folder in folders] + [('TOTAL', sum(counts.values()))])
            local_representatives = root / f'cleaned_{stem}.vsearch98.fasta'
            local_groups = root / 'vsearch_local_groups'
            local_summary = root / 'vsearch_local_groups.tsv'
            local_memberships = root / 'vsearch_local_memberships.tsv'
            with checkpoints.stage('vsearch_local', [local_representatives, local_summary, local_memberships],
                                   extra=[local_groups]) as execute:
                if execute:
                    cluster_local_groups(runner, local, local_representatives, local_groups,
                                         local_summary, local_memberships, args.threads)
            # References remain independent of local groups through clustering.
            reference_prefix = 'BOLD' if args.reference_source == 'bold' else 'COI'
            references = root / f'{reference_prefix}_{stem}.fasta'
            accessions = root / 'reference_accessions.txt'
            if args.reference_source == 'bold':
                download = root / 'references_bold.fasta'
                (root / 'reference_query.txt').write_text(f'BOLD header contains (case-insensitive): {name}\n')
                bold_metadata = root / 'bold_reference_metadata.tsv'
                with checkpoints.stage('bold_extraction', [download, accessions, bold_metadata]) as execute:
                    if execute:
                        count_bold = extract_bold(bold_files, name, download, bold_metadata, accessions)
                        runner.note(f'Extracted {count_bold} BOLD records matching {name}')
            else:
                query = reference_query(taxid, args.include_mitogenomes)
                (root / 'reference_query.txt').write_text(query + '\n')
                search = root / 'coi_search.xml'
                download = root / 'references_downloaded.fasta'
                with checkpoints.stage('download', [download, accessions], extra=[search], legacy=True) as execute:
                    if execute:
                        runner.run(['esearch', '-db', 'nuccore', '-query', query], search)
                        runner.run(['efetch', '-format', 'acc'], accessions, search)
                        if not accessions.read_text().strip():
                            raise ValueError('NCBI returned no COI accessions')
                        runner.run(['efetch', '-format', 'fasta'], download, search)
            with checkpoints.stage('reference_policy', [references, root / 'reference_filter.tsv', root / 'reference_accessions_retained.txt'], extra=(), legacy=True) as execute:
                if execute:
                    retained_refs, excluded_refs = filter_references(
                        download, references, root / 'reference_filter.tsv',
                        root / 'reference_accessions_retained.txt', args.include_mitogenomes)
                    runner.note(f'Reference policy: include_mitogenomes={args.include_mitogenomes}; '
                                f'retained {retained_refs}, excluded {excluded_refs} downloaded records')
            length_filtered = root / f'{reference_prefix}_{stem}.length500-2000.fasta'
            with checkpoints.stage('reference_lengths', [length_filtered], extra=(), legacy=True) as execute:
                if execute:
                    runner.run(['seqkit', 'seq', '--min-len', '500', '--max-len', '2000', references], length_filtered)
            reference_count = sum(1 for _ in fasta(length_filtered))
            if not reference_count:
                raise ValueError('No reference sequences remain in the inclusive 500–2000 nt range')
            runner.note(f'References retained after 500–2000 nt filtering: {reference_count}')
            reference_representatives = root / f'{reference_prefix}_{stem}.vsearch98.fasta'
            reference_uc = root / f'{reference_prefix}_{stem}.vsearch98.uc'
            with checkpoints.stage('vsearch_references', [reference_representatives, reference_uc],
                                   extra=[root / 'vsearch_references.stdout.log']) as execute:
                if execute:
                    vsearch_cluster(runner, length_filtered, reference_representatives, reference_uc,
                                    root / 'vsearch_references.stdout.log', args.threads)
            cluster_counts = [(dataset, sum(1 for _ in fasta(source)), sum(1 for _ in fasta(target)), '0.98')
                              for dataset, source, target in [('local', local, local_representatives),
                                                               ('reference', length_filtered, reference_representatives)]]
            tsv(root / 'clustering_counts.tsv', ['dataset', 'input_sequences', 'consensus_sequences', 'identity'], cluster_counts)
            # Combine consensus sets only after each source has been clustered.
            representatives = root / f'combined_{stem}.fasta'
            with checkpoints.stage('combined_representatives', [representatives]) as execute:
                if execute:
                    seen = set()
                    with representatives.open('x') as handle:
                        for source in (local_representatives, reference_representatives):
                            for header, sequence in fasta(source):
                                identifier = header.split()[0]
                                if identifier in seen or not sequence:
                                    raise ValueError(f'Duplicate ID or empty sequence: {identifier}')
                                seen.add(identifier)
                                write_record(handle, header, sequence)
                    runner.note(f'Combined {len(seen)} independently clustered representatives')
            # Nuclear/noncoding workflows can retain every consensus without translating it.
            if args.skip_translation:
                retained = root / f'combined_{slug}_{stage_label}_{args.date}.fasta'
                with checkpoints.stage('selection', [retained], legacy=False) as execute:
                    if execute:
                        shutil.copyfile(representatives, retained)
                runner.note('Translation and amino-acid selection skipped; retained all consensuses')
            else:
                translated = root / f'combined_{stem}.six_frames.faa'
                with checkpoints.stage('translation', [translated], extra=(), legacy=True) as execute:
                    if execute:
                        runner.run(['seqkit', 'translate', '-f', '6', '-F', '-T', args.translation_table,
                                    representatives], translated)
                retained = root / f'combined_{slug}_{stage_label}_{args.date}.fasta'
                with checkpoints.stage('retained', [retained, root / 'translation_lengths.tsv'],
                                       legacy_validator=lambda: validate_retained_checkpoint(
                                           representatives, translated, retained, root / 'translation_lengths.tsv',
                                           args.min_aa, args.aa_length_mode)) as execute:
                    if execute:
                        count = filter_translations(representatives, translated, retained,
                                                    root / 'translation_lengths.tsv', args.min_aa, args.aa_length_mode)
        count = sum(1 for _ in fasta(retained if retained.exists() else retained.with_suffix('.aln.fasta')))
        runner.note(f'Retained {count} nucleotide consensus sequences for alignment')
        if count < 4:
            raise ValueError('Fewer than four representatives passed; FASTA saved, phylogeny not run')
        alignment = retained.with_suffix('.aln.fasta')
        # kc-align publishes a completion marker only after mixed-mode alignment validates.
        if args.kcalign and not alignment_ready:
            from kcalign_alignment import align_kcalign
            kcalign_prefix = str(alignment).removesuffix('.aln.fasta')
            artifacts = [Path(kcalign_prefix+suffix) for suffix in (
                '.pre_align.aa.fasta', '.kcalign_selection.tsv', '.kcalign_ids.tsv', '.kcalign.stdout.log')]
            with checkpoints.stage('alignment', [alignment], extra=artifacts+[root/'alignment.complete.json'], legacy=False) as execute:
                if execute:
                    align_kcalign(runner, retained, alignment, args.kcalign_reference, args.translation_table, args.threads)
                    (root/'alignment.complete.json').write_text(json.dumps({
                        'algorithm': 'kc-align', 'alignment_sha256': file_digest(alignment),
                        'retained_sha256': file_digest(retained)}, indent=2)+'\n')
            alignment_ready = True
        if not alignment_ready:
            # A final-looking file without a completion record is an interrupted legacy result.
            if alignment.exists():
                runner.note(f'Removing interrupted alignment: {alignment}')
                alignment.unlink()
            with checkpoints.stage('alignment', [alignment, root / 'alignment_orientation.tsv'], extra=[root / 'alignment.complete.json'], legacy=False) as execute:
                if execute:
                    partial = Path(str(alignment) + '.partial')
                    partial.unlink(missing_ok=True)
                    try:
                        runner.run(['mafft', '--retree', '2', '--maxiterate', '0', '--adjustdirectionaccurately','--op', '20', '--ep', '18',
                                '--thread', args.threads, retained], partial)
                    except BaseException:
                        partial.unlink(missing_ok=True)
                        raise
                    # MAFFT prefixes reverse-complemented sequence IDs with _R_. Restore tracking headers.
                    aligned = list(fasta(partial))
                    originals = {header: sequence for header, sequence in fasta(retained)}
                    restored, orientation = [], []
                    for header, sequence in aligned:
                        reverse = header not in originals and header.startswith('_R_') and header[3:] in originals
                        original = header[3:] if reverse else header
                        if original not in originals:
                            raise ValueError(f'Unexpected MAFFT header: {header}')
                        restored.append((original, sequence))
                        orientation.append((original.split()[0], '-' if reverse else '+'))
                    if len(restored) != len(originals) or len({h for h, _ in restored}) != len(originals):
                        raise ValueError('MAFFT did not preserve all representative IDs')
                    if len({len(sequence) for _, sequence in restored}) != 1:
                        raise ValueError('MAFFT output is not an alignment')
                    with partial.open('w') as handle:
                        for header, sequence in restored:
                            write_record(handle, header, sequence)
                    tsv(root / 'alignment_orientation.tsv', ['sequence_id', 'alignment_strand'], orientation)
                    valid_alignment(partial, retained)
                    partial.replace(alignment)
                    (root / 'alignment.complete.json').write_text(json.dumps({
                        'algorithm': 'FFT-NS-2', 'alignment_sha256': file_digest(alignment),
                        'retained_sha256': file_digest(retained)}, indent=2) + '\n')
        prefix = root / f'combined_{slug}_{stage_label}_{args.date}'
        # Keep the untrimmed alignment so resume can start here without repeating alignment.
        pruned = Path(str(prefix) + '.clipkit.aln.fasta')
        trim_marker = root / 'clipkit.complete.json'
        trim_log = root / 'clipkit.stdout.log'
        trim_settings = {'mode': 'smart-gap', 'sequence_type': 'nt',
                         'input_sha256': file_digest(alignment)}
        if trim_marker.exists():
            try:
                if json.loads(trim_marker.read_text()) != trim_settings:
                    checkpoints.dirty = True
            except (ValueError, OSError):
                checkpoints.dirty = True
        with checkpoints.stage('clipkit', [pruned, trim_marker],
                               extra=[trim_log, Path(str(pruned)+'.partial')], legacy=False) as execute:
            if execute:
                # Persist invalidation before trimming, including runs created before this stage existed.
                checkpoints.invalidated.update(('phylogeny', 'figure'))
                checkpoints.invalidated_path.write_text(json.dumps(sorted(checkpoints.invalidated)))
                partial = Path(str(pruned)+'.partial')
                try:
                    runner.run(['clipkit', alignment, '-m', 'gappy', '-s', 'nt',
                                '-if', 'fasta', '-of', 'fasta', '-o', partial], trim_log)
                    validate_outputs([partial])
                    original_records = dict(fasta(alignment))
                    trimmed_records = list(fasta(partial))
                    if len(trimmed_records) != len(original_records) or {h for h, _ in trimmed_records} != set(original_records):
                        raise ValueError('ClipKIT changed sequence headers or the sequence set')
                    lengths = {len(seq) for _, seq in trimmed_records}
                    if len(lengths) != 1 or 0 in lengths or next(iter(lengths)) > len(next(iter(original_records.values()))):
                        raise ValueError('ClipKIT produced an empty or invalid alignment')
                    partial.replace(pruned)
                    trim_marker.write_text(json.dumps(trim_settings, indent=2)+'\n')
                finally:
                    partial.unlink(missing_ok=True)
        # Both alignment branches feed the same trimmed nucleotide phylogeny and figure stages.
        bootstrap_flag = '-b' if args.bootstrap_method == 'standard' else '-B'
        with checkpoints.stage('phylogeny', [Path(str(prefix) + '.treefile')], extra=[root / 'iqtree.stdout.log'] + [Path(str(prefix)+suffix) for suffix in ('.iqtree', '.log', '.ckp.gz', '.bionj', '.mldist', '.model.gz', '.contree', '.boottrees', '.ufboot', '.splits.nex')], legacy=False) as execute:
            if execute:
                runner.run(['iqtree', '-s', pruned, '-m', 'MFP', bootstrap_flag, '1000',
                            '-T', args.threads, '--seed', '666', '--prefix', prefix], root / 'iqtree.stdout.log')
                if not Path(str(prefix) + '.treefile').is_file():
                    raise ValueError('IQ-TREE did not produce the expected treefile')
        if args.reference_source == 'bold':
            reference_genbank = root / 'bold_reference_metadata.tsv'
            validate_outputs([reference_genbank])
        else:
            reference_genbank = root / 'reference_records.gb.xml'
            # Metadata is only needed for plotting, never to resume alignment or IQ-TREE.
            metadata_checkpoints = Checkpoints(runner, resume)
            with metadata_checkpoints.stage('reference_metadata', [reference_genbank], legacy=True) as execute:
                if execute:
                    runner.run(['efetch', '-db', 'nuccore', '-format', 'gb', '-mode', 'xml'],
                               reference_genbank, root / 'reference_accessions.txt')
            checkpoints.records['reference_metadata'] = metadata_checkpoints.records['reference_metadata']
        with checkpoints.stage('figure', [Path(str(prefix) + '.figure' + suffix) for suffix in ('.pdf', '.svg', '.midpoint.nwk', '.tips.tsv', '.plot.json')], extra=[root / 'tree_figure.stdout.log'], legacy=False) as execute:
            if execute:
                runner.run([sys.executable, Path(__file__).with_name('plot_phylogeny.py'),
                            '--tree', str(prefix) + '.treefile', '--analysis-dir', args.analysis_dir,
                            '--reference-metadata', reference_genbank,
                            '--output-prefix', str(prefix) + '.figure'], root / 'tree_figure.stdout.log')
        (root / 'postprocessing.complete').write_text('complete\n')
        (root / 'postprocessing.failed').unlink(missing_ok=True)
        runner.note(f'Completed: {prefix}.treefile')
    except Exception as error:
        for partial in root.glob('*.aln.fasta.partial'):
            partial.unlink(missing_ok=True)
        runner.note(f'FAILED: {error}')
        (root / 'postprocessing.failed').write_text(str(error) + '\n')
        raise


def main():
    """Parse source, filtering and alignment choices before starting or resuming a run."""
    parser = argparse.ArgumentParser(description=__doc__)
    reference_group = parser.add_mutually_exclusive_group(required=True)
    reference_group.add_argument('--download', dest='reference_source', action='store_const', const='download', help='Download NCBI COI references using Entrez Direct')
    reference_group.add_argument('--bold', dest='reference_source', action='store_const', const='bold', help='Select local BOLD FASTA headers containing the taxon name; no network')
    parser.add_argument('--bold-db', type=Path, help='BOLD gzipped FASTA file or directory (default: data/bold)')
    parser.add_argument('--include-mitogenomes', action='store_true',
                        help='Include matching whole mitogenome records (excluded by default)')
    parser.add_argument('--taxon', required=True, help='NCBI taxid or exact taxon name')
    parser.add_argument('--folders', nargs='+', required=True)
    parser.add_argument('--analysis-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True, help='New output directory, or existing run with --resume')
    parser.add_argument('--restart-alignment', action='store_true', help='With --resume, reuse retained sequences but rebuild alignment/tree/figure; permits changing aligner and trusted reference')
    parser.add_argument('--resume', action='store_true', help='Validate and reuse completed stages in the same output directory')
    parser.add_argument('--date', default=None, help='YYYYMMDD; defaults to today, or saved run date with --resume')
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--kcalign-reference', type=Path, help='Trusted single-sequence, in-frame coding FASTA required with --kcalign')
    parser.add_argument('--skip-translation', action='store_true', help='Skip six-frame translation and amino-acid selection; retain all consensuses (e.g. nuclear noncoding loci)')
    aligners = parser.add_mutually_exclusive_group()
    aligners.add_argument('--kcalign', action='store_true', help='Use kc-align mixed-mode codon alignment instead of MAFFT')
    aligners.add_argument('--mafft', dest='kcalign', action='store_false', help='Use MAFFT FFT-NS-2 (default)')
    parser.set_defaults(kcalign=False)
    parser.add_argument('--translation-table', type=int, default=5, help='SeqKit genetic code (5: invertebrate mitochondrial)')
    parser.add_argument('--min-aa', type=int, default=200)
    parser.add_argument('--aa-length-mode', choices=['stop-free'], default='stop-free',
                        help='Require consecutive amino acids without stops or unknown residues')
    parser.add_argument('--bootstrap-method', choices=['standard', 'ultrafast'], default='standard')
    args = parser.parse_args()
    if args.bold_db is not None and args.reference_source != 'bold':
        parser.error('--bold-db requires --bold')
    if min(args.threads, args.min_aa, args.translation_table) < 1:
        parser.error('Numeric settings must be positive')
    if args.date is not None and not re.fullmatch(r'\d{8}', args.date):
        parser.error('--date must have YYYYMMDD format')
    try:
        run_pipeline(args)
    except (ValueError, OSError, subprocess.CalledProcessError, ET.ParseError) as error:
        parser.exit(1, f'Post-processing failed: {error}\n')


if __name__ == '__main__':
    main()
