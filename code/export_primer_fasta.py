#!/usr/bin/env python3
"""Export cleaned reads of one barcode folder that carry a given forward tag (e.g. bc1087)."""
import argparse
import csv
import os
from pathlib import Path
import tempfile

from export_sintax_fasta import analysis_inputs, read_rows, write_selected
from sequence_headers import export_header

AUDIT = ['source_folder', 'readid', 'raw_ont_read_id', 'forward_PCR_primer', 'reverse_PCR_primer',
         'assigned_rank', 'assigned_label', 'assigned_confidence', 'taxonomy_status', 'contamination']


def export_primer(*, primer, barcode, analysis_dir, output, include_unclassified=False):
    output = Path(output)
    if output.suffix.lower() not in ('.fasta', '.fa', '.fna'):
        raise ValueError('--output must be an uncompressed .fasta, .fa or .fna file')
    destinations = [output, output.with_suffix('.reads.tsv'), output.with_suffix('.counts.tsv')]
    for path in destinations:
        if path.exists():
            raise FileExistsError(f'Refusing to overwrite {path}')
    folders = [Path(barcode).name] + (['unclassified'] if include_unclassified and barcode != 'unclassified' else [])
    output.parent.mkdir(parents=True, exist_ok=True)
    counts, total = [], 0
    with tempfile.TemporaryDirectory(prefix='.primer_export_', dir=output.parent) as directory:
        staged = [Path(directory)/p.name for p in destinations]
        with staged[0].open('w') as fasta_out, staged[1].open('w', newline='') as audit:
            writer = csv.writer(audit, delimiter='\t', lineterminator='\n')
            writer.writerow(AUDIT + ['export_header'])
            for folder in folders:
                table, clean = analysis_inputs(analysis_dir, folder)
                selected = {}
                for row in read_rows(table):
                    if row['forward_PCR_primer'] == primer:
                        header = export_header(folder, row['forward_PCR_primer'], row['reverse_PCR_primer'],
                                               row['assigned_label'], row['contamination'], row['readid'])
                        selected[row['readid']] = header
                        writer.writerow([row[c] for c in AUDIT] + [header])
                count = write_selected(fasta_out, clean, selected, folder)
                counts.append((folder, count))
                total += count
                print(f'{folder}: {count:,} cleaned reads with forward tag {primer}', flush=True)
        if not total:
            raise ValueError(f'No cleaned reads carry forward tag {primer} in {", ".join(folders)}')
        with staged[2].open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(['source_folder', 'reads_exported'])
            writer.writerows(counts + [('TOTAL', total)])
        published = []
        try:
            for stage, destination in zip(staged, destinations):
                os.link(stage, destination)
                published.append(destination)
        except BaseException:
            for path in published:
                path.unlink()
            raise
    print(f'Wrote {total:,} cleaned reads to {output}')
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--primer', required=True, help='Forward tag name, e.g. bc1087')
    parser.add_argument('--barcode', required=True, help='Sample folder, e.g. barcode15')
    parser.add_argument('--analysis-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--include-unclassified', action='store_true', help='Also search the unclassified folder')
    args = parser.parse_args()
    try:
        export_primer(**vars(args))
    except (ValueError, OSError) as error:
        parser.exit(1, f'Primer export failed: {error}\n')


if __name__ == '__main__':
    main()
