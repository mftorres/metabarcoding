"""Shared discovery for barcode folders and the ONT unclassified folder."""
import argparse
import json
from pathlib import Path
import re
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def project_path(path):
    """Store paths inside the project relative to it, so analyses survive moving the project."""
    path = Path(path).resolve()
    return str(path.relative_to(PROJECT_ROOT)) if path.is_relative_to(PROJECT_ROOT) else str(path)


def from_project(value):
    """Inverse of project_path."""
    path = Path(value)
    return path if path.is_absolute() else PROJECT_ROOT / path


def is_sample(name):
    return bool(re.fullmatch(r'barcode\d+', name)) or name == 'unclassified'


def dataset_path(value):
    path = Path(value).expanduser()
    if not path.exists() and len(path.parts) == 1:
        path = PROJECT_ROOT / 'data' / path
    if not path.is_dir():
        raise FileNotFoundError(f'Missing dataset or sample directory: {value}')
    return path.resolve()


def raw_files(folder):
    return sorted(p for p in Path(folder).iterdir()
                  if p.is_file() and p.name.endswith(('.fastq.gz', '.fq.gz')))


def cleaned_ready(folder):
    """Read processing output of the current pipeline: primer table plus cleaned reads."""
    return all((Path(folder)/'filtered'/name).is_file() for name in ('asinput.fastq.gz', 'primer_tags.tsv.gz'))


def discover_samples(values):
    result, names = [], set()
    for value in values:
        path = dataset_path(value)
        candidates = [path] if is_sample(path.name) else sorted(
            (p for p in path.iterdir() if p.is_dir() and is_sample(p.name)),
            key=lambda p: (p.name == 'unclassified', int(p.name[7:]) if p.name.startswith('barcode') else 0))
        for folder in candidates:
            if not (raw_files(folder) or cleaned_ready(folder)):
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
    root = from_project(params['merged_dir']).resolve()
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
    args = parser.parse_args()
    try:
        for folder in discover_samples(args.folders):
            print(folder)
    except (OSError, ValueError) as error:
        parser.exit(1, f'{error}\n')


if __name__ == '__main__':
    main()
