#!/usr/bin/env python3
"""Build full and reduced raw-read tracking tables, including an Insecta subset."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import shutil
import sys
from pathlib import Path


CONSENSUS_RE = re.compile(r"^consensus_(?:combinedrawreads|asinput)_(\d+_\d+)\(")
SUBGROUP_FASTA_RE = re.compile(r"^asinput_(\d+_\d+)\.fasta$")
SPLIT_READ_RE = re.compile(r"^split-by-adapter-(left|right)-(.+)$")


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


def read_amplicon_index(path: Path):
    internal_to_ont = {}
    processed_by_raw = {}
    processed_lengths = {}
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            internal = row["amplicon_sorter_read_number"]
            processed_id = row["ont_read_id"]
            internal_to_ont[internal] = processed_id
            processed_lengths[processed_id] = row["read_length"]

            split_match = SPLIT_READ_RE.match(processed_id)
            if split_match is None:
                fragment = "original"
                raw_id = processed_id
            else:
                fragment, raw_id = split_match.groups()

            processed_by_raw.setdefault(raw_id, []).append((processed_id, fragment))

    return internal_to_ont, processed_by_raw, processed_lengths


def read_subgroup_memberships(combined_dir: Path, internal_to_ont: dict[str, str]):
    memberships = {}
    duplicate_memberships = 0

    for path in sorted(combined_dir.glob("asinput_*_*.fasta")):
        match = SUBGROUP_FASTA_RE.match(path.name)
        if match is None:
            continue
        subgroup = match.group(1)
        for internal_id in fasta_ids(path):
            ont_id = internal_to_ont.get(internal_id)
            if ont_id is None:
                continue
            previous = memberships.setdefault(ont_id, subgroup)
            if previous != subgroup:
                duplicate_memberships += 1

    return memberships, duplicate_memberships


def read_kraken_assignments(kraken_path: Path, read_map_path: Path):
    kraken_taxids = {}
    with open(kraken_path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t", 4)
            if len(fields) >= 3:
                kraken_taxids[fields[1]] = fields[2]

    consensus_by_group = {}
    with open(read_map_path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            kraken_id = row["kraken_read_id"]
            group = row.get("amplicon_sorter_group", "")
            if not group:
                match = CONSENSUS_RE.match(kraken_id)
                if match is not None:
                    group = match.group(1)
            if group:
                consensus_by_group[group] = (
                    kraken_id,
                    kraken_taxids.get(kraken_id, ""),
                )

    return consensus_by_group


def read_target_kraken_membership(
    target_kraken_path: Path,
    read_map_path: Path,
    merge_folder: str,
    barcode: str,
):
    """Map target Kraken IDs back to consensus IDs."""
    target_kraken_ids = set()
    with open(target_kraken_path, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        required = {"readID", "merge_folder", "barcode"}
        missing = required.difference(header)
        if missing:
            raise ValueError(
                f"Missing columns in {target_kraken_path}: {sorted(missing)}. "
                "This table must be tab-delimited."
            )
        read_id_index = header.index("readID")
        merge_index = header.index("merge_folder")
        barcode_index = header.index("barcode")
        last_needed_index = max(read_id_index, merge_index, barcode_index)
        for line in handle:
            # Only split through the last needed field. The Kraken table can be
            # multi-gigabyte and its remaining taxonomy/primer columns are not
            # needed for this intersection.
            fields = line.rstrip("\r\n").split("\t", last_needed_index + 1)
            if len(fields) <= last_needed_index:
                continue
            if fields[merge_index] == merge_folder and fields[barcode_index] == barcode:
                target_kraken_ids.add(fields[read_id_index])

    target_consensus_ids = set()
    mapped_kraken_ids = set()
    with open(read_map_path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            kraken_id = row["kraken_read_id"]
            if kraken_id not in target_kraken_ids:
                continue
            mapped_kraken_ids.add(kraken_id)
            target_consensus_ids.add(kraken_id)

    unmapped_count = len(target_kraken_ids.difference(mapped_kraken_ids))
    return target_consensus_ids, len(target_kraken_ids), unmapped_count


def derived_output_paths(full_output: Path, barcode: str):
    if full_output.suffix != ".tsv":
        raise ValueError("--output must end in .tsv; the script creates the .tsv.gz file")

    reduced_output = full_output.with_name(f"{full_output.stem}_reduced.tsv")
    insecta_output = full_output.with_name(f"{barcode}_insecta_read_primer_group_taxid.tsv")
    insecta_reduced_output = insecta_output.with_name(f"{insecta_output.stem}_reduced.tsv")
    return reduced_output, insecta_output, insecta_reduced_output


def gzip_and_remove(path: Path):
    """Atomically create path.gz, then remove the uncompressed full table."""
    compressed = Path(f"{path}.gz")
    temporary = Path(f"{compressed}.tmp")
    with open(path, "rb") as source, gzip.open(
        temporary,
        "wb",
        compresslevel=6,
    ) as target:
        shutil.copyfileobj(source, target, length=1024 * 1024)
    temporary.replace(compressed)
    path.unlink()
    return compressed


def read_primer_hits(path: Path, merge_folder: str, barcode: str):
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
    return hits, fields


def find_single_kraken_output(combined_dir: Path) -> tuple[Path, Path]:
    kraken_outputs = sorted(combined_dir.glob("*consensussequences_noNUMT_*filter.fasta.kraken"))
    if len(kraken_outputs) != 1:
        raise FileNotFoundError(
            f"Expected one noNUMT Kraken output in {combined_dir}, found {len(kraken_outputs)}")

    kraken_output = kraken_outputs[0]
    kraken_input = Path(str(kraken_output)[:-len(".kraken")])
    kraken_map = kraken_input.with_name(f"{kraken_input.stem}.read_map.tsv")
    return kraken_output, kraken_map


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--barcode-dir", required=True, type=Path)
    parser.add_argument("--primer-best-hits", required=True, type=Path)
    parser.add_argument(
        "--insecta-kraken-reads",
        required=True,
        type=Path,
        help=(
            "Tab-delimited kraken_parser Insecta read table; its .csv extension "
            "may be retained"
        ),
    )
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    barcode_dir = args.barcode_dir.resolve()
    barcode = barcode_dir.name
    merge_folder = barcode_dir.parent.name
    combined_dir = barcode_dir / "filtered" / "combined"
    amplicon_map = combined_dir / "amplicon_sorter_read_map.tsv"
    kraken_output, kraken_map = find_single_kraken_output(combined_dir)

    for required in [
        amplicon_map,
        kraken_map,
        kraken_output,
        args.primer_best_hits,
        args.insecta_kraken_reads,
    ]:
        if not required.is_file():
            raise FileNotFoundError(required)

    print("Reading amplicon_sorter mappings...", file=sys.stderr)
    internal_to_ont, processed_by_raw, processed_lengths = read_amplicon_index(amplicon_map)
    memberships, duplicate_memberships = read_subgroup_memberships(combined_dir, internal_to_ont)
    print("Reading Kraken assignments...", file=sys.stderr)
    consensus_by_group = read_kraken_assignments(
        kraken_output,
        kraken_map,
    )
    print("Tracing Insecta Kraken IDs back to amplicon_sorter reads...", file=sys.stderr)
    (
        insecta_consensus_ids,
        insecta_kraken_id_count,
        unmapped_insecta_kraken_id_count,
    ) = read_target_kraken_membership(
        args.insecta_kraken_reads,
        kraken_map,
        merge_folder,
        barcode,
    )
    print("Reading barcode primer best hits...", file=sys.stderr)
    primer_hits, primer_fields = read_primer_hits(args.primer_best_hits, merge_folder, barcode)

    columns = [
        "merge_folder",
        "barcode",
        "raw_ont_read_id",
        "raw_fastq_file",
        "raw_read_length",
        *primer_fields,
        "entered_amplicon_sorter",
        "processed_read_id",
        "processed_fragment",
        "processed_read_length",
        "amplicon_sorter_group",
        "amplicon_sorter_group_consensus_id",
        "group_consensus_kraken_taxid",
    ]

    reduced_columns = [
        "raw_ont_read_id",
        "forward_patternName",
        "reverse_patternName",
        "group_consensus_kraken_taxid",
    ]
    reduced_output, insecta_output, insecta_reduced_output = derived_output_paths(
        args.output,
        barcode,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    raw_count = 0
    full_row_count = 0
    insecta_row_count = 0
    insecta_raw_ids = set()
    with (
        open(args.output, "w", encoding="utf-8", newline="") as full_handle,
        open(reduced_output, "w", encoding="utf-8", newline="") as reduced_handle,
        open(insecta_output, "w", encoding="utf-8", newline="") as insecta_handle,
        open(insecta_reduced_output, "w", encoding="utf-8", newline="") as insecta_reduced_handle,
    ):
        writer_options = {
            "delimiter": "\t",
            "extrasaction": "ignore",
            "lineterminator": "\n",
        }
        full_writer = csv.DictWriter(full_handle, fieldnames=columns, **writer_options)
        reduced_writer = csv.DictWriter(
            reduced_handle,
            fieldnames=reduced_columns,
            **writer_options,
        )
        insecta_writer = csv.DictWriter(insecta_handle, fieldnames=columns, **writer_options)
        insecta_reduced_writer = csv.DictWriter(
            insecta_reduced_handle,
            fieldnames=reduced_columns,
            **writer_options,
        )
        for writer in [full_writer, reduced_writer, insecta_writer, insecta_reduced_writer]:
            writer.writeheader()

        for raw_id, fastq_name, read_length in raw_fastq_records(barcode_dir):
            raw_count += 1
            processed_records = processed_by_raw.get(raw_id) or [("", "")]
            for processed_id, fragment in processed_records:
                group = memberships.get(processed_id, "")
                consensus_id, consensus_taxid = consensus_by_group.get(group, ("", ""))
                row = {
                    "merge_folder": merge_folder,
                    "barcode": barcode,
                    "raw_ont_read_id": raw_id,
                    "raw_fastq_file": fastq_name,
                    "raw_read_length": read_length,
                    "entered_amplicon_sorter": "yes" if processed_id else "no",
                    "processed_read_id": processed_id,
                    "processed_fragment": fragment,
                    "processed_read_length": processed_lengths.get(processed_id, ""),
                    "amplicon_sorter_group": group,
                    "amplicon_sorter_group_consensus_id": consensus_id,
                    "group_consensus_kraken_taxid": consensus_taxid,
                }
                row.update(primer_hits.get(raw_id, {}))
                full_writer.writerow(row)
                reduced_writer.writerow(row)
                full_row_count += 1

                is_insecta = consensus_id in insecta_consensus_ids
                if is_insecta:
                    insecta_writer.writerow(row)
                    insecta_reduced_writer.writerow(row)
                    insecta_row_count += 1
                    insecta_raw_ids.add(raw_id)

    print(f"Wrote {full_row_count:,} rows from {raw_count:,} raw reads to {args.output}", file=sys.stderr)
    print(
        f"Wrote {insecta_row_count:,} Insecta rows from "
        f"{len(insecta_raw_ids):,} raw reads to {insecta_output}",
        file=sys.stderr,
    )
    print(f"Wrote reduced tables to {reduced_output} and {insecta_reduced_output}", file=sys.stderr)
    print(
        f"Insecta Kraken IDs for {merge_folder}/{barcode}: {insecta_kraken_id_count:,}; "
        f"unmapped in the Kraken read map: {unmapped_insecta_kraken_id_count:,}",
        file=sys.stderr,
    )
    print(
        f"Processed read/fragment records linked to raw reads: "
        f"{sum(len(records) for records in processed_by_raw.values()):,}",
        file=sys.stderr,
    )
    print(f"Primer best-hit rows loaded for this barcode: {len(primer_hits):,}", file=sys.stderr)
    print(f"Reads assigned to an amplicon_sorter subgroup: {len(memberships):,}", file=sys.stderr)
    if duplicate_memberships:
        print(
            f"Warning: {duplicate_memberships:,} reads occurred in more than one subgroup; first kept",
            file=sys.stderr,
        )

    print("Compressing full tracking tables...", file=sys.stderr)
    compressed_full = gzip_and_remove(args.output)
    compressed_insecta = gzip_and_remove(insecta_output)
    print(f"Wrote {compressed_full}", file=sys.stderr)
    print(f"Wrote {compressed_insecta}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
