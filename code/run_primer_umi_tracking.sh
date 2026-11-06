#!/usr/bin/env bash
set -euo pipefail

CONFIG_FILE="${CONFIG_FILE:-config/pipeline.env}"

if [[ -f "$CONFIG_FILE" ]]; then
  # shellcheck disable=SC1090
  source "$CONFIG_FILE"
fi

: "${PRIMER_UMI_FASTA:=data/saproxdna_primertags_seq5_2026_1.fasta}"
: "${PRIMER_UMI_MAX_MISMATCH:=10}"
: "${PRIMER_UMI_SEQKIT_THREADS:=8}"
: "${MAX_JOBS:=2}"
: "${FORCE:=0}"

usage() {
  cat >&2 <<'EOF'
Usage:
  bash code/run_primer_umi_tracking.sh FOLDER [FOLDER ...]

Each FOLDER can be either:
  - a parent folder containing barcodeNN and unclassified folders, for example data/sequencing2
  - a single barcodeNN or unclassified folder, for example data/sequencing2/barcode01

Examples:
  bash code/run_primer_umi_tracking.sh data/sequencing2
  MAX_JOBS=4 bash code/run_primer_umi_tracking.sh data/sequencing2 data/sequencing3
  bash code/run_primer_umi_tracking.sh data/sequencing2/barcode01

Useful variables:
  PRIMER_UMI_FASTA=data/saproxdna_primertags_seq5_2026_1.fasta
  PRIMER_UMI_MAX_MISMATCH=5
  PRIMER_UMI_SEQKIT_THREADS=8
  MAX_JOBS=2
  FORCE=1
EOF
}

require_command() {
  local command_name="$1"
  if ! command -v "$command_name" >/dev/null 2>&1; then
    echo "Missing command: $command_name" >&2
    exit 1
  fi
}

run_barcode() {
  local barcode_dir="$1"
  local barcode_name
  local raw_read_list
  local output
  local partial_output
  local log
  local file_count

  barcode_name="$(basename "$barcode_dir")"
  raw_read_list="$barcode_dir/rawreads_${barcode_name}.txt"
  output="$barcode_dir/${barcode_name}_primer_hits.tsv"
  partial_output="${output}.partial.$$"
  log="$barcode_dir/${barcode_name}_primer_hits.log"

  echo "[$barcode_dir] listing raw FASTQ files"
  find "$barcode_dir" -maxdepth 1 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' \) | sort -V > "$raw_read_list"
  file_count="$(wc -l < "$raw_read_list")"

  if [[ "$file_count" -eq 0 ]]; then
    echo "[$barcode_dir] no raw *fastq.gz files found" >&2
    return 1
  fi

  if [[ "$FORCE" != "1" && -s "$output" ]]; then
    echo "[$barcode_dir] primer UMI hits already exist: $output"
    return 0
  fi

  rm -f "$partial_output"
  echo "[$barcode_dir] seqkit locate on $file_count FASTQ files"
  {
    printf 'Command:'
    printf ' %q' \
      seqkit locate \
      --max-mismatch "$PRIMER_UMI_MAX_MISMATCH" \
      --threads "$PRIMER_UMI_SEQKIT_THREADS" \
      --ignore-case \
      --immediate-output \
      --pattern-file "$PRIMER_UMI_FASTA" \
      --infile-list "$raw_read_list"
    printf ' > %q\n\n' "$partial_output"
  } > "$log"
  # seems like seqkit locate doesnt have an option to limit the numnber of hits per subject-query or per query
  # relaxed, large max-mismatches values will produce huge files that need to be loaded by chuncks for the notebook
  # also, seqkit deals with positive and negative strand search by default, e.g. reverse complements of the primers+tags
  seqkit locate \
    --max-mismatch "$PRIMER_UMI_MAX_MISMATCH" \
    --threads "$PRIMER_UMI_SEQKIT_THREADS" \
    --ignore-case \
    --immediate-output \
    --pattern-file "$PRIMER_UMI_FASTA" \
    --infile-list "$raw_read_list" \
    > "$partial_output" 2>> "$log"

  [[ -s "$partial_output" ]] || {
    echo "[$barcode_dir] seqkit produced an empty output: $partial_output" | tee -a "$log" >&2
    rm -f "$partial_output"
    return 1
  }

  mv "$partial_output" "$output"
  echo "[$barcode_dir] wrote $output"
}

main() {
  if [[ "$#" -lt 1 ]]; then
    usage
    exit 2
  fi

  [[ "$MAX_JOBS" =~ ^[1-9][0-9]*$ ]] || { echo "MAX_JOBS must be positive" >&2; exit 2; }
  require_command seqkit

  [[ -f "$PRIMER_UMI_FASTA" ]] || {
    echo "Missing primer UMI FASTA: $PRIMER_UMI_FASTA" >&2
    exit 1
  }

  local barcode_dir discovered
  discovered="$(python3 "$(dirname "${BASH_SOURCE[0]}")/sample_folders.py" --require-raw "$@")" || exit 1
  local barcode_folders=()
  mapfile -t barcode_folders <<< "$discovered"

  local failed=0 pid
  local pids=()
  for barcode_dir in "${barcode_folders[@]}"; do
    run_barcode "$barcode_dir" &
    pids+=("$!")
    if [[ "${#pids[@]}" -ge "$MAX_JOBS" ]]; then
      for pid in "${pids[@]}"; do wait "$pid" || failed=1; done
      pids=()
    fi
  done
  for pid in "${pids[@]}"; do wait "$pid" || failed=1; done

  exit "$failed"
}

main "$@"
