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
  - a parent folder containing barcode* folders, for example data/mergedseq5
  - a single barcode folder, for example data/mergedseq5/barcode01

Examples:
  bash code/run_primer_umi_tracking.sh data/mergedseq5
  MAX_JOBS=4 bash code/run_primer_umi_tracking.sh data/mergedseq5 data/mergedseq6
  bash code/run_primer_umi_tracking.sh data/mergedseq5/barcode01

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

discover_barcode_folders() {
  local input="$1"

  [[ -d "$input" ]] || {
    echo "Missing folder: $input" >&2
    return 1
  }

  if [[ "$(basename "$input")" == barcode* ]]; then
    printf '%s\n' "$input"
  else
    find "$input" -mindepth 1 -maxdepth 1 -type d -name 'barcode*' | sort -V
  fi
}

run_barcode() {
  local barcode_dir="$1"
  local barcode_name
  local raw_read_list
  local output
  local log
  local file_count

  barcode_name="$(basename "$barcode_dir")"
  raw_read_list="$barcode_dir/rawreads_${barcode_name}.txt"
  output="$barcode_dir/${barcode_name}_primer_hits.tsv"
  log="$barcode_dir/${barcode_name}_primer_hits.log"

  echo "[$barcode_dir] listing raw FASTQ files"
  find "$barcode_dir" -maxdepth 1 -type f -name '*fastq.gz' | sort -V > "$raw_read_list"
  file_count="$(wc -l < "$raw_read_list")"

  if [[ "$file_count" -eq 0 ]]; then
    echo "[$barcode_dir] no raw *fastq.gz files found" >&2
    return 1
  fi

  if [[ "$FORCE" != "1" && -s "$output" ]]; then
    echo "[$barcode_dir] primer UMI hits already exist: $output"
    return 0
  fi

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
    printf ' > %q\n\n' "$output"
  } > "$log"

  seqkit locate \
    --max-mismatch "$PRIMER_UMI_MAX_MISMATCH" \
    --threads "$PRIMER_UMI_SEQKIT_THREADS" \
    --ignore-case \
    --immediate-output \
    --pattern-file "$PRIMER_UMI_FASTA" \
    --infile-list "$raw_read_list" \
    > "$output" 2>> "$log"

  [[ -s "$output" ]] || {
    echo "[$barcode_dir] seqkit produced an empty output: $output" | tee -a "$log" >&2
    return 1
  }

  echo "[$barcode_dir] wrote $output"
}

wait_for_slot() {
  while [[ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]]; do
    sleep 2
  done
}

main() {
  if [[ "$#" -lt 1 ]]; then
    usage
    exit 2
  fi

  require_command seqkit

  [[ -f "$PRIMER_UMI_FASTA" ]] || {
    echo "Missing primer UMI FASTA: $PRIMER_UMI_FASTA" >&2
    exit 1
  }

  local input
  local barcode_dir
  local barcode_folders=()

  for input in "$@"; do
    while IFS= read -r barcode_dir; do
      barcode_folders+=("$barcode_dir")
    done < <(discover_barcode_folders "$input")
  done

  if [[ "${#barcode_folders[@]}" -eq 0 ]]; then
    echo "No barcode* folders found" >&2
    exit 1
  fi

  local failed=0
  for barcode_dir in "${barcode_folders[@]}"; do
    wait_for_slot
    run_barcode "$barcode_dir" &
  done

  for job in $(jobs -rp); do
    if ! wait "$job"; then
      failed=1
    fi
  done

  exit "$failed"
}

main "$@"
