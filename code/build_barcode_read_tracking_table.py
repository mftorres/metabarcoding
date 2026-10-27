#!/usr/bin/env python3
"""Join original ONT reads, primer hits and direct SINTAX assignments."""
import argparse
import csv
import gzip
from pathlib import Path
from contextlib import ExitStack
from sintax_parser import COLUMNS, parse_output
from sintax_analysis import primer_name_map

def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def fasta_ids(path: Path):
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                yield line[1:].split()[0]


def raw_fastq_records(barcode_dir: Path):
    paths = sorted(barcode_dir.glob("*.fastq.gz")) + sorted(barcode_dir.glob("*.fq.gz"))
    if not paths:
        raise FileNotFoundError(f"No raw *.fastq.gz or *.fq.gz files found in {barcode_dir}")

    for path in paths:
        with open_text(path) as handle:
            record_number = 0
            while True:
                header = handle.readline()
                if not header:
                    break
                sequence = handle.readline()
                plus = handle.readline()
                quality = handle.readline()
                record_number += 1
                if not sequence or not plus or not quality:
                    raise ValueError(f"Incomplete FASTQ record {record_number} in {path}")
                if not header.startswith("@"):
                    raise ValueError(f"Invalid FASTQ header at record {record_number} in {path}")
                yield header[1:].split()[0], path.name, len(sequence.rstrip("\r\n"))


def read_primer_hits(path: Path, merge_folder: str, barcode: str, mapping):
    fields = [
        "forward_patternName",
        "forward_strand",
        "forward_start",
        "forward_end",
        "forward_mismatch_count",
        "forward_primer_identity",
        "reverse_patternName",
        "reverse_strand",
        "reverse_start",
        "reverse_end",
        "reverse_mismatch_count",
        "reverse_primer_identity",
    ]
    hits = {}
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["merged_folder"] == merge_folder and row["barcode_folder"] == barcode:
                hits[row["seqID"]] = {field: row.get(field, "") for field in fields}
    canonical = set(mapping.values())
    for hit in hits.values():
        for direction in ('forward', 'reverse'):
            original = hit.pop(direction + '_patternName', '')
            if original and original not in mapping and original not in canonical:
                raise ValueError(f'Unmapped primer name: {original}')
            hit[direction + '_PCR_primer'] = mapping.get(original, original)
    fields = [field.replace('_patternName', '_PCR_primer') for field in fields]
    return hits, fields



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--barcode-dir', required=True, type=Path)
    parser.add_argument('--primer-best-hits', required=True, type=Path)
    parser.add_argument('--sintax-output', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    reference_root = Path(__file__).resolve().parents[1] / 'data'
    parser.add_argument('--mft-primers', type=Path, default=reference_root/'saproxdna_primertags_seq5_2026_1.fasta')
    parser.add_argument('--pcr-primers', type=Path, default=reference_root/'PCRsaproxdna_primertags_seq5_2026_1.fasta')
    args = parser.parse_args()
    folder = args.barcode_dir.resolve()
    barcode, merge_folder = folder.name, folder.parent.name
    if args.sintax_output.parent.name != barcode:
        raise ValueError('SINTAX output folder must match the sample folder name')
    mapping = primer_name_map(args.mft_primers, args.pcr_primers)
    hits, primer_fields = read_primer_hits(args.primer_best_hits, merge_folder, barcode, mapping)
    assignments = {}
    seen = set()
    for row in parse_output(args.sintax_output, cutoff=0.8):
        if row['readid'] in seen:
            raise ValueError(f"Duplicate SINTAX query ID: {row['readid']}")
        seen.add(row['readid'])
        assignments.setdefault(row['raw_ont_read_id'], []).append(row)
    if args.output.suffix != '.tsv':
        raise ValueError('--output must end in .tsv')
    reduced = args.output.with_name(args.output.stem + '_reduced.tsv')
    insecta = args.output.with_name(barcode + '_insecta_read_sintax_reduced.tsv')
    columns = ['merge_folder', 'raw_fastq_file', 'raw_read_length', *COLUMNS, *primer_fields]
    reduced_columns = ['merge_folder', 'source_folder', 'barcode', 'source_file', 'source_name',
                       'raw_ont_read_id', 'readid', 'forward_PCR_primer',
                       'reverse_PCR_primer', 'sintax_taxid', 'sintax_name',
                       'taxlabel', 'taxlabel_lowest', 'taxlabel_lowestrank', 'accepted_class',
                       'accepted_family', 'accepted_genus', 'accepted_species']
    for path in (args.output, reduced, insecta):
        if path.exists():
            raise FileExistsError(f'Refusing to overwrite {path}')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with ExitStack() as stack:
        writers = []
        for path, fields in ((args.output, columns), (reduced, reduced_columns), (insecta, reduced_columns)):
            handle = stack.enter_context(path.open('x', encoding='utf-8', newline=''))
            writer = csv.DictWriter(handle, fields, delimiter='\t', extrasaction='ignore', lineterminator='\n')
            writer.writeheader()
            writers.append(writer)
        for raw_id, filename, length in raw_fastq_records(folder):
            for assignment in assignments.get(raw_id, [{}]):
                row = dict(assignment)
                row.update(merge_folder=merge_folder, raw_fastq_file=filename,
                           raw_read_length=length, raw_ont_read_id=raw_id,
                           source_folder=barcode, barcode=barcode,
                           source_file=str(args.sintax_output.resolve()), source_name=args.sintax_output.name)
                row.update(hits.get(raw_id, {}))
                writers[0].writerow(row)
                writers[1].writerow(row)
                if row.get('accepted_class') == 'Insecta_50557':
                    writers[2].writerow(row)
                count += 1
    print(f'Wrote {count:,} read/fragment rows to {args.output}')

if __name__ == '__main__':
    main()
