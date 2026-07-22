#!/usr/bin/env python3
"""Build a table linking Kraken consensus record IDs to amplicon_sorter groups."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


CONSENSUS_RE = re.compile(r"^consensus_(?:combinedrawreads|asinput)_(\d+_\d+)\(")


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


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Map Kraken consensus record IDs back to amplicon_sorter groups."
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--targets", required=True, nargs="+", type=Path)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)

    columns = [
        "merge_folder",
        "barcode",
        "kraken_read_id",
        "amplicon_sorter_group",
        "kraken_sequence_length",
    ]

    with open(args.output, "w", encoding="utf-8") as out:
        out.write("\t".join(columns) + "\n")
        for target in args.targets:
            merge_folder, barcode = barcode_metadata(target)
            group_name = target.stem

            for kraken_read_id, sequence in fasta_records(target):
                match = CONSENSUS_RE.match(kraken_read_id)
                group_name = match.group(1) if match is not None else target.stem

                row = [
                    merge_folder,
                    barcode,
                    kraken_read_id,
                    group_name,
                    str(len(sequence)),
                ]
                out.write("\t".join(row) + "\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
