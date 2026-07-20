#!/usr/bin/env bash
set -euo pipefail

CONFIG_FILE="${CONFIG_FILE:-config/pipeline.env}"

if [[ -f "$CONFIG_FILE" ]]; then
  # shellcheck disable=SC1090
  source "$CONFIG_FILE"
fi

: "${NONUMT_MIN_AA_LENGTH:=100}"
: "${NONUMT_TRANSLATION_TABLE:=5}"
: "${NONUMT_SEQKIT_THREADS:=6}"
: "${MAX_JOBS:=2}"
: "${FORCE:=0}"

usage() {
  cat >&2 <<'EOF'
Usage:
  bash code/filter_consensus_noNUMT_seqkit.sh FOLDER [FOLDER ...]

Each FOLDER can be either:
  - a parent folder containing barcode* folders, for example data/mergedseq5
  - a single barcode folder, for example data/mergedseq5/barcode01

Examples:
  bash code/filter_consensus_noNUMT_seqkit.sh data/mergedseq5
  NONUMT_MIN_AA_LENGTH=200 bash code/filter_consensus_noNUMT_seqkit.sh data/mergedseq5
  FORCE=1 bash code/filter_consensus_noNUMT_seqkit.sh data/mergedseq5/barcode01

Useful variables:
  CONFIG_FILE=config/pipeline.env
  NONUMT_MIN_AA_LENGTH=100
  NONUMT_TRANSLATION_TABLE=5
  NONUMT_SEQKIT_THREADS=6
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
  local combined_dir="$barcode_dir/filtered/combined"
  local consensus_fasta
  local input_name
  local input_stem
  local output_fasta
  local input_sequence_count
  local output_sequence_count
  local output_count=0

  [[ -d "$combined_dir" ]] || {
    echo "[$barcode_dir] missing combined folder: $combined_dir" >&2
    return 1
  }

  shopt -s nullglob
  local consensus_fastas=("$combined_dir"/*consensussequences.fasta)
  shopt -u nullglob

  if [[ "${#consensus_fastas[@]}" -eq 0 ]]; then
    echo "[$barcode_dir] no *consensussequences.fasta files found" >&2
    return 0
  fi

  for consensus_fasta in "${consensus_fastas[@]}"; do
    [[ "$consensus_fasta" != *noNUMT* ]] || continue

    input_name="$(basename "$consensus_fasta")"
    input_stem="${input_name%.fasta}"
    output_fasta="$combined_dir/${input_stem}_aa_noNUMT_${NONUMT_MIN_AA_LENGTH}filter.fasta"
    input_sequence_count="$(grep -c '^>' "$consensus_fasta" || true)"

    if [[ "$FORCE" != "1" && -s "$output_fasta" ]]; then
      echo "[$barcode_dir] noNUMT output already exists: $output_fasta"
      output_sequence_count="$(grep -c '^>' "$output_fasta" || true)"
      echo "[$barcode_dir] sequence counts: input=$input_sequence_count output=$output_sequence_count"
      output_count=$((output_count + 1))
      continue
    fi

    echo "[$barcode_dir] translating and filtering $consensus_fasta"
    echo "[$barcode_dir] writing $output_fasta"

    seqkit translate "$consensus_fasta" \
      --frame 6 \
      --append-frame \
      -T "$NONUMT_TRANSLATION_TABLE" \
      -j "$NONUMT_SEQKIT_THREADS" \
      --trim \
      | seqkit seq \
        --min-len "$NONUMT_MIN_AA_LENGTH" \
        -g \
        -o "$output_fasta"

    output_sequence_count="$(grep -c '^>' "$output_fasta" || true)"
    echo "[$barcode_dir] sequence counts: input=$input_sequence_count output=$output_sequence_count"

    output_count=$((output_count + 1))
  done

  echo "[$barcode_dir] noNUMT amino-acid FASTA files written or already present: $output_count"
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
