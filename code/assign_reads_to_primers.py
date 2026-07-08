#!/usr/bin/env python3
import argparse
import csv
import gzip
from pathlib import Path


IUPAC = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "U": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("GC"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}

COMPLEMENT = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn", "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")


def read_fasta(path):
    records = []
    name = None
    seq_parts = []
    with open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq_parts).upper().replace("U", "T")))
                name = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if name is not None:
        records.append((name, "".join(seq_parts).upper().replace("U", "T")))
    return records


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1].upper().replace("U", "T")


def compile_primer(name, seq, strand):
    allowed = [IUPAC.get(base, IUPAC["N"]) for base in seq]
    informative = sum(1 for base in seq if base != "N")
    return {
        "name": name,
        "seq": seq,
        "strand": strand,
        "allowed": allowed,
        "length": len(seq),
        "informative": informative or len(seq),
    }


def build_seed_index(primers, seed_length):
    seed_index = {}
    if seed_length <= 0:
        return seed_index
    for primer_index, primer in enumerate(primers):
        seq = primer["seq"]
        for offset in range(0, len(seq) - seed_length + 1):
            seed = seq[offset : offset + seed_length]
            if all(base in "ACGT" for base in seed):
                seed_index.setdefault(seed, []).append((primer_index, offset))
    return seed_index


def iter_fastq_gz(path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq = handle.readline().strip().upper().replace("U", "T")
            _plus = handle.readline()
            _qual = handle.readline()
            if not _qual:
                return
            yield header[1:].strip().split()[0], seq


def candidate_starts(read_len, primer_len, scan_bases):
    max_start = read_len - primer_len
    if max_start < 0:
        return []
    if scan_bases <= 0 or scan_bases >= read_len:
        return range(max_start + 1)

    starts = set()
    left_end = min(scan_bases, read_len)
    for start in range(max(0, left_end - primer_len + 1)):
        starts.add(start)

    right_start = max(0, read_len - scan_bases)
    for start in range(right_start, max_start + 1):
        starts.add(start)

    return sorted(starts)


def score_at(read, primer, start):
    matches = 0
    comparable = 0
    for offset, allowed in enumerate(primer["allowed"]):
        base = read[start + offset]
        if base not in IUPAC:
            continue
        if primer["seq"][offset] != "N":
            comparable += 1
        if base in allowed:
            matches += 1
    denominator = primer["informative"]
    return matches, comparable, matches / denominator if denominator else 0.0


def best_primer_match(read, primers, scan_bases, seed_index, seed_length):
    best = None
    candidates = set()

    if seed_index:
        if scan_bases <= 0 or scan_bases >= len(read):
            spans = [(0, len(read))]
        else:
            spans = [(0, min(scan_bases, len(read))), (max(0, len(read) - scan_bases), len(read))]

        for span_start, span_end in spans:
            for read_offset in range(span_start, max(span_start, span_end - seed_length + 1)):
                seed = read[read_offset : read_offset + seed_length]
                for primer_index, primer_offset in seed_index.get(seed, []):
                    start = read_offset - primer_offset
                    primer = primers[primer_index]
                    if 0 <= start <= len(read) - primer["length"]:
                        candidates.add((primer_index, start))
    else:
        for primer_index, primer in enumerate(primers):
            for start in candidate_starts(len(read), primer["length"], scan_bases):
                candidates.add((primer_index, start))

    for primer_index, start in candidates:
        primer = primers[primer_index]
        matches, comparable, score = score_at(read, primer, start)
        mismatches = primer["informative"] - matches
        hit = (
            score,
            matches,
            -mismatches,
            primer["name"],
            primer["strand"],
            start,
            start + primer["length"],
            primer["length"],
            comparable,
            mismatches,
        )
        if best is None or hit > best:
            best = hit
    return best


def main():
    parser = argparse.ArgumentParser(
        description="Assign raw FASTQ reads to the best matching primer+UMI sequence."
    )
    parser.add_argument("--fastq-dir", required=True, type=Path, help="Folder containing *.fastq.gz files.")
    parser.add_argument("--primers", required=True, type=Path, help="Primer+UMI FASTA file.")
    parser.add_argument("--output", required=True, type=Path, help="Output TSV path.")
    parser.add_argument(
        "--scan-bases",
        type=int,
        default=150,
        help="Scan this many bases from each read end. Use 0 to scan full reads. Default: 150.",
    )
    parser.add_argument(
        "--min-score",
        type=float,
        default=0.0,
        help="Minimum best score to report a primer name. Lower scores are reported as unassigned.",
    )
    parser.add_argument(
        "--seed-length",
        type=int,
        default=10,
        help="Exact A/C/G/T seed length used to find candidate primer placements. Use 0 for exhaustive search.",
    )
    args = parser.parse_args()

    primer_records = read_fasta(args.primers)
    primers = []
    for name, seq in primer_records:
        primers.append(compile_primer(name, seq, "+"))
        primers.append(compile_primer(name, reverse_complement(seq), "-"))
    seed_index = build_seed_index(primers, args.seed_length)

    fastq_files = sorted(args.fastq_dir.glob("*.fastq.gz"))
    if not fastq_files:
        raise FileNotFoundError(f"No *.fastq.gz files found in {args.fastq_dir}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "wt", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(
            [
                "read_id",
                "best_primer_umi",
                "score",
                "matches",
                "mismatches",
                "strand",
                "start_0based",
                "end_0based",
                "primer_length",
                "fastq_file",
            ]
        )
        for fastq in fastq_files:
            for read_id, read in iter_fastq_gz(fastq):
                best = best_primer_match(read, primers, args.scan_bases, seed_index, args.seed_length)
                if best is None or best[0] < args.min_score:
                    writer.writerow([read_id, "unassigned", "", "", "", "", "", "", "", fastq.name])
                    continue
                score, matches, _neg_mismatches, name, strand, start, end, length, _comparable, mismatches = best
                writer.writerow(
                    [
                        read_id,
                        name,
                        f"{score:.6f}",
                        matches,
                        mismatches,
                        strand,
                        start,
                        end,
                        length,
                        fastq.name,
                    ]
                )


if __name__ == "__main__":
    main()
