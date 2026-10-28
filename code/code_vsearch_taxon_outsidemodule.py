# conda activate taxonomy
# just run this on dry in the command line
python - \
  data/exports/Ceratopogonidae_barcodes01-21_seq5.fasta \
  data/exports/Ceratopogonidae_centroids98 <<'PY'
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, "code")
from postprocess_sintax import fasta, write_record, tsv

source = Path(sys.argv[1])
output_dir = Path(sys.argv[2])
threads = 2

# Read only the already-extracted taxon dataset.
groups = defaultdict(list)
for header, sequence in fasta(source):
    fields = header.split()[0].split("-", 5)
    if len(fields) != 6:
        raise ValueError(f"Unexpected tracking header: {header}")
    groups[tuple(fields[:3])].append((header, sequence))

if not groups:
    raise ValueError("Input FASTA contains no sequences")

# Refuse to overwrite an existing output directory.
output_dir.mkdir(parents=True, exist_ok=False)
combined = output_dir / "Ceratopogonidae.centroids98.fasta"
summary = []

with combined.open("x") as combined_handle:
    for number, (group, records) in enumerate(sorted(groups.items()), 1):
        prefix = output_dir / f"group{number:06d}"
        group_input = prefix.with_suffix(".input.fasta")
        centroids = prefix.with_suffix(".centroids.fasta")
        membership = prefix.with_suffix(".uc")

        with group_input.open("x") as handle:
            for header, sequence in records:
                write_record(handle, header, sequence)

        # Same identity/strand settings as Module 2; request centroids.
        command = [
            "vsearch", "--cluster_fast", str(group_input),
            "--id", "0.98", "--iddef", "2",
            "--strand", "both", "--qmask", "none",
            "--minseqlength", "1",
            "--maxseqlength", str(max(len(seq) for _, seq in records)),
            "--notrunclabels",
            "--centroids", str(centroids),
            "--uc", str(membership),
            "--threads", str(threads)]

        with prefix.with_suffix(".log").open("x") as log:
            subprocess.run(command, stdout=log, stderr=log, check=True)

        count = 0
        for header, sequence in fasta(centroids):
            write_record(combined_handle, header, sequence)
            count += 1

        summary.append((prefix.name, *group, len(records), count))

tsv(output_dir / "clustering_summary.tsv",
    ["group_id", "barcode", "forward_primer", "reverse_primer",
     "input_reads", "centroids"], summary)

print(f"Combined centroids: {combined}")
print(f"Input reads: {sum(row[-2] for row in summary)}")
print(f"Centroids: {sum(row[-1] for row in summary)}")
PY