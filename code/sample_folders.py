"""Shared discovery for barcode folders and the ONT unclassified folder."""
import argparse
import json
from pathlib import Path
import re
import sys


def is_sample(name):
    return bool(re.fullmatch(r'barcode\d+', name)) or name == 'unclassified'


def dataset_path(value):
    path = Path(value).expanduser()
    if not path.exists() and len(path.parts) == 1:
        path = Path(__file__).resolve().parents[1] / 'data' / path
    if not path.is_dir():
        raise FileNotFoundError(f'Missing dataset or sample directory: {value}')
    return path.resolve()


def raw_files(folder):
    return sorted(p for p in Path(folder).iterdir()
                  if p.is_file() and p.name.endswith(('.fastq.gz', '.fq.gz')))


def discover_samples(values, require_raw=False):
    result, names = [], set()
    for value in values:
        path = dataset_path(value)
        candidates = [path] if is_sample(path.name) else sorted(
            (p for p in path.iterdir() if p.is_dir() and is_sample(p.name)),
            key=lambda p: (p.name == 'unclassified', int(p.name[7:]) if p.name.startswith('barcode') else 0))
        for folder in candidates:
            usable = bool(raw_files(folder)) or (not require_raw and (folder/'filtered/asinput.fastq.gz').is_file())
            if not usable:
                print(f'Skipping {folder}: no direct raw FASTQ or usable cleaned input; merge pass/fail first if needed', file=sys.stderr)
                continue
            if folder.name in names:
                raise ValueError(f'Duplicate sample name: {folder.name}; use separate runs for separate datasets')
            result.append(folder.resolve()); names.add(folder.name)
    if not result:
        raise ValueError('No nonempty barcodeNN or unclassified sample folders found')
    return result


def analysis_folders(values, analysis_dir):
    """Expand a parent using completed analysis metadata, including unclassified."""
    params = json.loads((Path(analysis_dir)/'analysis_parameters.json').read_text())
    available = params['analysed_folders']
    root = Path(params['merged_dir']).resolve()
    result = []
    for value in values or available:
        if value in available:
            selected = [value]
        else:
            if is_sample(str(value)):
                raise ValueError(f'{value} was not analysed; regenerate Read processing to include it')
            path = dataset_path(value)
            if path == root:
                selected = [name for name in available if (root/name).is_dir()]
            elif path.parent == root and path.name in available:
                selected = [path.name]
            else:
                raise ValueError(f'{value} is not part of the completed analysis at {analysis_dir}')
        for name in selected:
            if name in result:
                raise ValueError(f'Duplicate sample: {name}')
            result.append(name)
    if not result:
        raise ValueError('No analysed folders selected')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folders', nargs='+')
    parser.add_argument('--require-raw', action='store_true')
    args = parser.parse_args()
    try:
        for folder in discover_samples(args.folders, args.require_raw):
            print(folder)
    except (OSError, ValueError) as error:
        parser.exit(1, f'{error}\n')


if __name__ == '__main__':
    main()
