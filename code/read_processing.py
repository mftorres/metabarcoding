#!/usr/bin/env python3
"""Module 1, Read processing: primers/tags, cleaning, SINTAX and report tables for one dataset."""
import argparse
from datetime import datetime
import os
from pathlib import Path
import subprocess
import sys

from sample_folders import dataset_path, discover_samples, is_sample


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('dataset', help='Parent folder of the barcode folders, e.g. sequencing2 or data/sequencing2')
    parser.add_argument('--output-dir', type=Path, help='New analysis directory; default: data/read_processing/DATASET_TIMESTAMP')
    parser.add_argument('--sintax-root', type=Path, help='Classification outputs; default: data/sintax/DATASET')
    parser.add_argument('--dry-run', action='store_true', help='Show what would run, without running it')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    dataset = dataset_path(args.dataset)
    if is_sample(dataset.name):
        parser.error('Provide the parent folder that contains the barcode folders')
    folders = discover_samples([dataset])
    sintax_root = (args.sintax_root or root/'data/sintax'/dataset.name).resolve()
    output = (args.output_dir or root/'data/read_processing'/f'{dataset.name}_{datetime.now():%Y%m%d_%H%M%S}').resolve()
    if output.exists():
        parser.error(f'Use a new output directory: {output}')
    command = ['bash', str(root/'code/run_metabarcoding_pipeline.sh')] + (['--dry-run'] if args.dry_run else []) + [str(dataset)]
    env = dict(os.environ, SINTAX_OUTPUT_ROOT=str(sintax_root), PYTHON=sys.executable)
    subprocess.run(command, cwd=root, env=env, check=True)
    print(f'Analysis of {", ".join(p.name for p in folders)} -> {output}/analysis', flush=True)
    if args.dry_run:
        return
    from sintax_analysis import run_analysis
    run_analysis(merged_dir=dataset, sintax_root=sintax_root, output_dir=output/'analysis',
                 folders=[p.name for p in folders], allow_missing=True)
    print(f'Read processing complete. Postprocessing --analysis-dir: {output}/analysis')


if __name__ == '__main__':
    main()
