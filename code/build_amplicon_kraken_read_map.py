#!/usr/bin/env python3
"""Build a table linking amplicon_sorter/Kraken record IDs to ONT read IDs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_amplicon_sorter_map(path: Path) -> dict[str, tuple[str, int]]:
    read_map: dict[str, tuple[str, int]] = {}
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"amplicon_sorter_read_number", "ont_read_id", "read_length"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(sorted(missing))}")

        for row in reader:
            read_map[row["amplicon_sorter_read_number"]] = (
                row["ont_read_id"],
                int(row["read_length"]),
            )

    return read_map


def open_text(path: Path):
    return open(path, "rt", encoding="utf-8")


def fasta_records(path: Path):
    record_id = None
    sequence_parts: list[str] = []

    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if record_id is not None:
                    yield record_id, "".join(sequence_parts)
                record_id = line[1:].split()[0]
                sequence_parts = []
            else:
                sequence_parts.append(line)

    if record_id is not None:
        yield record_id, "".join(sequence_parts)


def barcode_metadata(path: Path) -> tuple[str, str]:
    parts = path.parts
    barcode_index = next((i for i, part in enumerate(parts) if part.startswith("barcode")), None)
    barcode = parts[barcode_index] if barcode_index is not None else ""
    merge_folder = parts[barcode_index - 1] if barcode_index and barcode_index > 0 else ""
    return merge_folder, barcode


def record_type(record_id: str) -> str:
    if record_id.isdigit():
        return "read"
    if record_id.startswith("consensus"):
        return "consensus"
    return "other"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Map amplicon_sorter FASTA/Kraken record IDs back to ONT read IDs."
    )
    parser.add_argument("--amplicon-sorter-read-map", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--targets", required=True, nargs="+", type=Path)
    args = parser.parse_args()

    read_map = read_amplicon_sorter_map(args.amplicon_sorter_read_map)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    columns = [
        "merge_folder",
        "barcode",
        "kraken_read_id",
        "record_type",
        "amplicon_sorter_file",
        "amplicon_sorter_group",
        "amplicon_sorter_read_number",
        "ont_read_id",
        "ont_read_length",
        "kraken_sequence_length",
    ]

    with open(args.output, "w", encoding="utf-8") as out:
        out.write("\t".join(columns) + "\n")
        for target in args.targets:
            merge_folder, barcode = barcode_metadata(target)
            group_name = target.stem

            for kraken_read_id, sequence in fasta_records(target):
                rec_type = record_type(kraken_read_id)
                amplicon_read_number = kraken_read_id if rec_type == "read" else ""
                ont_read_id = ""
                ont_read_length = ""

                if amplicon_read_number:
                    original = read_map.get(amplicon_read_number)
                    if original is not None:
                        ont_read_id, read_length = original
                        ont_read_length = str(read_length)

                row = [
                    merge_folder,
                    barcode,
                    kraken_read_id,
                    rec_type,
                    target.name,
                    group_name,
                    amplicon_read_number,
                    ont_read_id,
                    ont_read_length,
                    str(len(sequence)),
                ]
                out.write("\t".join(row) + "\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
