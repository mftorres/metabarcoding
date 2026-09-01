conda activate saproxdna

# check paths throughout the metabarcoding folder
# merge failed and passed reads by barcode
# create a folder ../data/merged/
# create folders "barcode[0-9]+" from 1 to 24 or 96
# mv *fastq.gz from passed/failed to corresponding ../data/merged/barcode
# set -euo pipefail
# mkdir -p data/merged

# for i in $(seq 1 96); do
#   mkdir -p "data/merged/$(printf 'barcode%02d' "$i")"
# done

moved=0

while IFS= read -r -d '' file; do
  source_barcode=$(basename "$(dirname "$file")")
  barcode_num="${source_barcode#barcode}"
  dest_barcode="$(printf 'barcode%02d' "$barcode_num")"

  dest="data/merged/${dest_barcode}/$(basename "$file")"

  if [ -e "$dest" ]; then
    printf 'collision: %s -> %s\n' "$file" "$dest" >&2
    exit 1
  fi

  cp "$file" "../data/merged/${dest_barcode}/"
  moved=$((moved + 1))
done < <(find data \( -path '*/fastq_pass/barcode*/*.fastq.gz' -o -path '*/fastq_fail/barcode*/*.fastq.gz' \) -print0)

printf 'moved %s files\n' "$moved"

# I could pass a fasta file with adapter and see if I can force fastplong to
# output sequences trimmed by adaptor sequences
# running from metabarcoding folder
FOLDER="data/merged/barcode59"

mkdir -p "$FOLDER"/filtered
mkdir -p "$FOLDER"/filtered/logs

# todo: generate file-based html reports
# ONT adapter sequences are not public
# --disable_adapter_trimmings to detect and trimm
# should work independently of passing a file to -a
# add a file to redirect stdout to for checking errors and adapters identified
# add check for empty folder from *gz or *fasta
for file in "$FOLDER"/*.fastq.gz; do
  [ -e "$file" ] || continue
  filename="$(basename "$file")"
  base="${filename%.fastq.gz}"
  out="${base}_filtered.fastq.gz"
  log="$FOLDER/filtered/logs/${base}.log"
  fastplong -i "$file" -o "$FOLDER/filtered/$out" \
      --unqualified_percent_limit 50 \
      --qualified_quality_phred 10 \
      --mean_qual 6 \
      --length_required 150 \
      --adapter_fasta data/ONT_24_native_barcode_sequences.fasta \
      -h "$FOLDER/filtered/${base}_report.html" \
      > "$log" 2>&1
done

# conda install bioconda::python-edlib
# conda install biopython
# conda install matplotlib

# get all filtered reads in a single file
cat "$FOLDER"/filtered/*_filtered.fastq.gz > "$FOLDER"/filtered/asinput.fastq.gz

# Metagenetic sample with several amplicons 
# unknown number of species, 30000 reads in the inputfile
# one low abundant species (< 2% reads):
# -maxr, --maxreads: Maximum number of reads to process. Default=10000
# -ar, --allreads: Use all reads from the inputfile between length limits. 
# This argument is still limited with --maxreads to have a hard limit for large files.
# -sg, --similar_genes: Similarity to sort genes in groups 
# (value between 50 and 100). Default=80.0
# -ssg, --similar_species_groups: Similarity to CREATE species groups (value between 50 and 100). Default=Estimate
# -ss, --similar_species: Similarity to ADD sequences to a species group (value between 50 and 100). Default=85.0
# -sc, --similar_consensus: Similarity to COMBINE groups based on the consensus sequence (value between 50 and 100). Default=96.0

python3 code/amplicon_sorter_mod.py -i "$FOLDER"/filtered/asinput.fastq.gz \
  -o "$FOLDER"/filtered/ \
  -min 300 -max 100000000 -ar -maxr 100000000 -np 12 \
  -ra -sg 99.0 -ss 99.0 -sc 99.0

# Prefix consensus IDs with the barcode so IDs remain unique when barcode
# consensus FASTAs are combined for downstream classification.
CONSENSUS_FASTA="$FOLDER/filtered/combined/asinput_consensussequences.fasta"
if [[ -s "$CONSENSUS_FASTA" ]]; then
  tmp_consensus="${CONSENSUS_FASTA}.prefixing.$$"
  barcode_name="$(basename "$FOLDER")"
  awk -v prefix="${barcode_name}_" \
    '/^>/ && substr($0, 2, length(prefix)) != prefix { sub(/^>/, ">" prefix) } { print }' \
    "$CONSENSUS_FASTA" > "$tmp_consensus"
  mv "$tmp_consensus" "$CONSENSUS_FASTA"
fi

# to kraken2, we only pass consensus sequences:
# asinput_consensussequences.fasta

# search for consensus files
# asinput_consensussequences.fasta
## collect matching files into an array for processing with kraken2
TARGET_FILES=()
for pattern in "$FOLDER/filtered/combined/"*consensussequences.fasta; do
  for f in $pattern; do
    [ -e "$f" ] || continue
    [[ "$f" != *noNUMT* ]] || continue
    TARGET_FILES+=("$f")
  done
done
# Example: print files and (optionally) run kraken2 via Docker on each
for file in "${TARGET_FILES[@]}"; do
  echo "Found: $file"
  hostdir=$(realpath "$(dirname "$file")")
  docker run --rm -v "$hostdir":/data \
  -v /mnt/d/kraken2DB/k2_core_nt_2025:/db staphb/kraken2 kraken2 \
  --db /db \
  --memory-mapping \
  --confidence 0.8 \
  --threads 10 \
  /data/"$(basename "$file")" \
  --output /data/"$(basename "$file").kraken" \
  --report /data/"$(basename "$file").report"
done
