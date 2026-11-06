#!/usr/bin/env python3
"""Combine ONT pass/fail gzip FASTQs per barcode/unclassified, preserving originals."""
import argparse
import json
import os
from pathlib import Path
import shutil
import tempfile
from sample_folders import dataset_path, is_sample

CATEGORIES = ('pass', 'fail', 'passed', 'failed', 'fastq_pass', 'fastq_fail')


def merge_dataset(dataset, dry_run=False):
    root = dataset_path(dataset)
    sources = {}
    # Support sample/pass/*.fastq.gz and pass/sample/*.fastq.gz layouts.
    for child in sorted(root.iterdir()):
        if not child.is_dir():
            continue
        if is_sample(child.name):
            for category in CATEGORIES:
                base = child/category
                if base.is_dir():
                    sources.setdefault(child.name, []).extend(p for p in base.rglob('*') if p.is_file() and p.name.endswith(('.fastq.gz', '.fq.gz')))
        elif child.name in CATEGORIES:
            for sample in sorted(child.iterdir()):
                if sample.is_dir() and is_sample(sample.name):
                    sources.setdefault(sample.name, []).extend(p for p in sample.rglob('*') if p.is_file() and p.name.endswith(('.fastq.gz', '.fq.gz')))
    plans = []
    for name, paths in sorted(sources.items()):
        paths = sorted(set(p.resolve() for p in paths))
        if not paths:
            continue
        folder = root/name
        output = folder/'pass_fail_merged.fastq.gz'
        manifest = folder/'pass_fail_merge.json'
        fingerprints = [{'path': str(p), 'bytes': p.stat().st_size, 'mtime_ns': p.stat().st_mtime_ns} for p in paths]
        total = sum(p.stat().st_size for p in paths)
        direct = [p for p in folder.glob('*') if p.is_file() and p.name.endswith(('.fastq.gz', '.fq.gz'))]
        if any(p != output for p in direct):
            raise ValueError(f'{folder}: direct raw FASTQ files already exist; refusing to risk counting the same reads twice')
        if output.exists():
            if manifest.is_file() and json.loads(manifest.read_text()) == {'sources': fingerprints, 'output_bytes': total} and output.stat().st_size == total:
                print(f'Reuse {output}: unchanged sources')
                continue
            raise ValueError(f'{output} already exists or source files changed; use a fresh dataset directory')
        if (folder/'filtered').exists() or (folder/f'{name}_primer_hits.tsv').exists():
            raise ValueError(f'{folder}: analyses already exist; merge into a fresh dataset before analysis')
        plans.append((folder, output, manifest, paths, fingerprints, total))
    if not sources:
        print('No nested pass/fail FASTQs found; existing flat sample folders need no merge.')
    for folder, output, manifest, paths, fingerprints, total in plans:
        print(f'{folder.name}: combine {len(paths)} files -> {output}')
        if dry_run:
            continue
        folder.mkdir(parents=True, exist_ok=True)
        # Concatenated gzip members remain a valid gzip stream; no decompression needed.
        with tempfile.NamedTemporaryFile(dir=folder, prefix='.merge-', delete=False) as handle:
            temporary = Path(handle.name)
            try:
                for source in paths:
                    with source.open('rb') as incoming:
                        shutil.copyfileobj(incoming, handle)
                handle.flush()
                if temporary.stat().st_size != total:
                    raise ValueError('Unexpected merge size')
                os.link(temporary, output)  # Exclusive publication; never overwrite a raw input.
            finally:
                temporary.unlink(missing_ok=True)
        with manifest.open('x') as handle:
            json.dump({'sources': fingerprints, 'output_bytes': total}, handle, indent=2)
            handle.write('\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('dataset', help='Dataset parent, e.g. sequencing2 or data/sequencing2')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    try:
        merge_dataset(args.dataset, args.dry_run)
    except (OSError, ValueError) as error:
        parser.exit(1, f'Merge failed: {error}\n')


if __name__ == '__main__':
    main()
