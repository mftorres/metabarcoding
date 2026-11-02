#!/usr/bin/env python3
"""Read processing: cleaning, SINTAX, primer tracking and reports for a dataset parent."""
import argparse
from datetime import datetime
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile

from sample_folders import dataset_path, discover_samples, is_sample


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('dataset', help='Parent folder, e.g. sequencing2 or data/sequencing2')
    parser.add_argument('--output-dir', type=Path, help='New run directory; default: data/read_processing/DATASET_TIMESTAMP')
    parser.add_argument('--sintax-root', type=Path, help='Classification outputs; default: data/sintax/DATASET (separate datasets)')
    parser.add_argument('--plate-map', type=Path, default=Path('saproxDNA/species_info_seq5/seq5_pcr_to_pools.tsv'))
    parser.add_argument('--dry-run', action='store_true', help='Print the steps without running or writing analyses')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    dataset = dataset_path(args.dataset)
    if is_sample(dataset.name):
        parser.error('Provide the parent containing the sample folders')
    folders = discover_samples([dataset])
    args.sintax_root = args.sintax_root or root/'data/sintax'/dataset.name
    output = (args.output_dir or root/'data/read_processing'/f'{dataset.name}_{datetime.now():%Y%m%d_%H%M%S_%f}').resolve()
    if output.exists():
        parser.error(f'Use a new output directory: {output}')
    commands = [
        ['bash', str(root/'code/run_metabarcoding_pipeline.sh'), str(dataset)],
        ['bash', str(root/'code/run_primer_umi_tracking.sh'), str(dataset)],
    ]
    for command in commands:
        print(shlex.join(command), flush=True)
    print(f'Parse primer hits for {", ".join(p.name for p in folders)} into {output}/primertrack', flush=True)
    print(f'Run post-classification notebooks into {output}/analysis', flush=True)
    if args.dry_run:
        return
    import nbformat
    from nbclient import NotebookClient
    output.mkdir(parents=True)
    env = dict(os.environ, SINTAX_OUTPUT_ROOT=str(args.sintax_root.resolve()))
    for command in commands:
        subprocess.run(command, cwd=root, env=env, check=True)
    # Use the launching environment for the primer-parser notebook's kernel.
    with tempfile.TemporaryDirectory(prefix='read_processing_kernel_') as temp:
        temp = Path(temp)
        spec = temp/'kernels/read_processing'; spec.mkdir(parents=True)
        (spec/'kernel.json').write_text(json.dumps({'argv': [sys.executable, '-m', 'ipykernel_launcher', '-f', '{connection_file}'],
                                                   'display_name': 'Read processing', 'language': 'python'}))
        os.environ.update(JUPYTER_PATH=str(temp)+os.pathsep+os.environ.get('JUPYTER_PATH',''),
                          MPLBACKEND='Agg', MPLCONFIGDIR=str(temp/'matplotlib'),
                          SEQKIT_INPUT_FOLDER=str(dataset), PRIMER_OUTPUT_DIR=str(output/'primertrack'))
        notebook = nbformat.read(root/'code/primer_umi_seqkit_parser.ipynb', as_version=4)
        try:
            NotebookClient(notebook, timeout=-1, kernel_name='read_processing', resources={'metadata': {'path': str(root)}}).execute()
        finally:
            nbformat.write(notebook, output/'primer_tracking.ipynb')
    best = list((output/'primertrack').glob('*_seqkit_primer_umi_best_hits_*.txt'))
    if len(best) != 1:
        raise ValueError('Expected one prepared best-hit table')
    subprocess.run([sys.executable, str(root/'code/run_sintax_notebooks.py'), str(dataset),
                    '--output-dir', str(output/'analysis'), '--executed-dir', str(output/'executed_notebooks'), '--primer-hits', str(best[0]),
                    '--sintax-root', str(args.sintax_root.resolve()), '--plate-map', str(args.plate_map.resolve())],
                   cwd=root, env=env, check=True)
    print(f'Read processing complete. Postprocessing --analysis-dir: {output}/analysis')


if __name__ == '__main__':
    main()
