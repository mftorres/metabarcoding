#!/usr/bin/env python3
"""Annotate noNUMT consensus FASTAs with Kraken names and trust status."""

from __future__ import annotations

import argparse
import csv
import os
import re
import tempfile
from pathlib import Path


DEFAULT_NAMES_DMP = Path("/mnt/d/kraken2DB/k2_core_nt_2025/names.dmp")
TRUE_VALUES = {"1", "true", "yes", "y"}


def sanitize_header(value: str) -> str:
    """Replace every non-alphanumeric run with one underscore."""
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_") or "unnamed"


def parse_taxid(value: str) -> int | None:
    value = str(value).strip()
    if not value:
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def find_single(barcode_dir: Path, pattern: str, description: str) -> Path:
    matches = sorted(barcode_dir.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(
            f"{barcode_dir.name}: expected one {description} matching {pattern!r}, "
            f"found {len(matches)}: {[str(path) for path in matches]}"
        )
    return matches[0]


def read_kraken_assignments(path: Path) -> tuple[dict[str, int | None], set[int]]:
    assignments: dict[str, int | None] = {}
    classified_taxids: set[int] = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected at least 3 tab-separated fields")
            status, sequence_id, taxid_text = fields[:3]
            taxid = parse_taxid(taxid_text)
            if status != "C" or taxid in {None, 0}:
                taxid = None
            else:
                classified_taxids.add(taxid)
            if sequence_id in assignments:
                raise ValueError(f"{path}:{line_number}: duplicate Kraken sequence ID {sequence_id!r}")
            assignments[sequence_id] = taxid
    return assignments, classified_taxids


def read_trusted_taxids(path: Path) -> set[int]:
    trusted_taxids: set[int] = set()
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"taxid_detected", "trusted"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing required columns {sorted(missing)}")
        for row in reader:
            if str(row["trusted"]).strip().lower() not in TRUE_VALUES:
                continue
            taxid = parse_taxid(row["taxid_detected"])
            if taxid not in {None, 0}:
                trusted_taxids.add(taxid)
    return trusted_taxids


def read_scientific_names(path: Path, wanted_taxids: set[int]) -> dict[int, str]:
    names: dict[int, str] = {}
    if not wanted_taxids:
        return names
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = [field.strip() for field in line.split("|")]
            if len(fields) < 4 or fields[3] != "scientific name":
                continue
            taxid = parse_taxid(fields[0])
            if taxid in wanted_taxids:
                names[taxid] = fields[1]
                if len(names) == len(wanted_taxids):
                    break
    return names


def fasta_records(path: Path):
    header: str | None = None
    sequence_lines: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if header is not None:
                    yield header, sequence_lines
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"{path}:{line_number}: empty FASTA header")
                sequence_lines = []
            elif header is None:
                if line:
                    raise ValueError(f"{path}:{line_number}: sequence before first FASTA header")
            else:
                sequence_lines.append(line)
    if header is not None:
        yield header, sequence_lines


def annotate_barcode(
    barcode_dir: Path,
    fasta_path: Path,
    kraken_path: Path,
    trusted_path: Path,
    taxon_names: dict[int, str],
    output_path: Path,
) -> dict[str, int]:
    assignments, _ = read_kraken_assignments(kraken_path)
    trusted_taxids = read_trusted_taxids(trusted_path)
    counts = {"total": 0, "trusted": 0, "notrusted": 0, "noid": 0}
    seen_headers: set[str] = set()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{output_path.name}.", suffix=".partial", dir=output_path.parent
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as output:
            for original_header, sequence_lines in fasta_records(fasta_path):
                sequence_id = original_header.split(maxsplit=1)[0]
                taxid = assignments.get(sequence_id)
                if taxid is None:
                    suffix = "noid"
                    counts["noid"] += 1
                else:
                    scientific_name = taxon_names.get(taxid)
                    if not scientific_name:
                        raise ValueError(
                            f"{barcode_dir.name}: scientific name missing for Kraken taxID {taxid}"
                        )
                    trust_status = "trusted" if taxid in trusted_taxids else "notrusted"
                    suffix = f"{sanitize_header(scientific_name)}_{trust_status}"
                    counts[trust_status] += 1

                new_header = sanitize_header(
                    f"{original_header}_{barcode_dir.name}_{suffix}"
                )
                if new_header in seen_headers:
                    raise ValueError(
                        f"{barcode_dir.name}: duplicate sanitized FASTA header {new_header!r}"
                    )
                seen_headers.add(new_header)
                output.write(f">{new_header}\n")
                for sequence_line in sequence_lines:
                    output.write(f"{sequence_line}\n")
                counts["total"] += 1
        os.replace(temporary_name, output_path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise
    return counts


def concatenate_fastas(input_paths: list[Path], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{output_path.name}.", suffix=".partial", dir=output_path.parent
    )
    seen_headers: set[str] = set()
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as output:
            for input_path in input_paths:
                with input_path.open(encoding="utf-8") as source:
                    for line in source:
                        if line.startswith(">"):
                            header = line[1:].strip()
                            if header in seen_headers:
                                raise ValueError(
                                    f"duplicate header in combined FASTA: {header!r}"
                                )
                            seen_headers.add(header)
                        output.write(line)
        os.replace(temporary_name, output_path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--merged-dir", type=Path, default=Path("data/mergedseq5"))
    parser.add_argument("--barcode-start", type=int, default=1)
    parser.add_argument("--barcode-end", type=int, default=24)
    parser.add_argument("--names-dmp", type=Path, default=DEFAULT_NAMES_DMP)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    barcode_dirs = [
        args.merged_dir / f"barcode{number:02d}"
        for number in range(args.barcode_start, args.barcode_end + 1)
    ]
    for barcode_dir in barcode_dirs:
        if not barcode_dir.is_dir():
            raise FileNotFoundError(f"missing barcode directory: {barcode_dir}")
    if not args.names_dmp.is_file():
        raise FileNotFoundError(f"missing NCBI names file: {args.names_dmp}")

    inputs = []
    all_classified_taxids: set[int] = set()
    for barcode_dir in barcode_dirs:
        if (barcode_dir / "filtered" / "no_consensus.done").is_file():
            inputs.append((barcode_dir, None, None, None))
            continue
        combined_dir = barcode_dir / "filtered" / "combined"
        fasta_path = find_single(
            combined_dir, "*consensussequences_noNUMT_*filter.fasta", "noNUMT consensus FASTA"
        )
        kraken_path = Path(f"{fasta_path}.kraken")
        if not kraken_path.is_file():
            raise FileNotFoundError(f"{barcode_dir.name}: missing Kraken output {kraken_path}")
        trusted_path = barcode_dir / f"{barcode_dir.name}_target_species_detection_trusted.tsv"
        if not trusted_path.is_file():
            raise FileNotFoundError(f"{barcode_dir.name}: missing trusted table {trusted_path}")
        _, classified_taxids = read_kraken_assignments(kraken_path)
        all_classified_taxids.update(classified_taxids)
        inputs.append((barcode_dir, fasta_path, kraken_path, trusted_path))

    taxon_names = read_scientific_names(args.names_dmp, all_classified_taxids)
    missing_names = sorted(all_classified_taxids - set(taxon_names))
    if missing_names:
        raise ValueError(
            f"{len(missing_names)} classified Kraken taxIDs have no scientific name in "
            f"{args.names_dmp}: {missing_names[:20]}"
        )

    barcode_outputs: list[Path] = []
    grand_counts = {"total": 0, "trusted": 0, "notrusted": 0, "noid": 0}
    for barcode_dir, fasta_path, kraken_path, trusted_path in inputs:
        output_path = (
            barcode_dir / "phylo" / f"{barcode_dir.name}_consensus_noNUMT_phylo.fasta"
        )
        if fasta_path is None:
            concatenate_fastas([], output_path)
            counts = {"total": 0, "trusted": 0, "notrusted": 0, "noid": 0}
        else:
            counts = annotate_barcode(
                barcode_dir,
                fasta_path,
                kraken_path,
                trusted_path,
                taxon_names,
                output_path,
            )
        barcode_outputs.append(output_path)
        for key in grand_counts:
            grand_counts[key] += counts[key]
        print(
            f"{barcode_dir.name}: total={counts['total']:,} "
            f"trusted={counts['trusted']:,} notrusted={counts['notrusted']:,} "
            f"noid={counts['noid']:,} -> {output_path}"
        )

    combined_output = (
        args.merged_dir
        / "consensus"
        / f"{args.merged_dir.name}_consensus_noNUMT_phylo.fasta"
    )
    concatenate_fastas(barcode_outputs, combined_output)
    print(
        f"Combined: total={grand_counts['total']:,} "
        f"trusted={grand_counts['trusted']:,} "
        f"notrusted={grand_counts['notrusted']:,} "
        f"noid={grand_counts['noid']:,} -> {combined_output}"
    )


if __name__ == "__main__":
    main()
