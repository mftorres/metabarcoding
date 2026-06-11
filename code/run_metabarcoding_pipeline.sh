#!/usr/bin/env bash
set -euo pipefail

CONFIG_FILE="${CONFIG_FILE:-config/pipeline.env}"

if [[ -f "$CONFIG_FILE" ]]; then
  # shellcheck disable=SC1090
  source "$CONFIG_FILE"
fi

: "${ADAPTER_FASTA:=data/ONT_24_native_barcode_sequences.fasta}"
: "${FASTPLONG_UNQUALIFIED_PERCENT_LIMIT:=50}"
: "${FASTPLONG_QUALIFIED_QUALITY_PHRED:=6}"
: "${FASTPLONG_MEAN_QUAL:=6}"
: "${FASTPLONG_LENGTH_REQUIRED:=150}"
: "${AMPLICON_SORTER_SCRIPT:=code/amplicon_sorter.py}"
: "${AMPLICON_THREADS:=14}"
: "${AMPLICON_MINLENGTH:=300}"
: "${AMPLICON_MAXLENGTH:=100000000}"
: "${AMPLICON_MAXREADS:=10000000000000000}"
: "${AMPLICON_SG:=99.0}"
: "${AMPLICON_SS:=99.0}"
: "${AMPLICON_SC:=99.0}"
: "${AMPLICON_EXTRA_ARGS:=}"
: "${KRAKEN_THREADS:=4}"
: "${KRAKEN_DB:=/mnt/d/kraken2DB/k2_core_nt_2025}"
: "${KRAKEN_CONFIDENCE:=0}"
: "${KRAKEN_MEMORY_MAPPING:=1}"
: "${KRAKEN_USE_DOCKER:=0}"
: "${KRAKEN_DOCKER_IMAGE:=staphb/kraken2}"
: "${MAX_JOBS:=2}"
: "${FORCE:=0}"
: "${FORCE_KRAKEN:=0}"

usage() {
  cat >&2 <<'EOF'
Usage:
  bash code/run_metabarcoding_pipeline.sh FOLDER [FOLDER ...]

Examples:
  bash code/run_metabarcoding_pipeline.sh data/merged/barcode58
  MAX_JOBS=3 bash code/run_metabarcoding_pipeline.sh data/merged/barcode58 data/merged/barcode60 data/merged/barcode61

Optional config:
  CONFIG_FILE=config/pipeline.env

Useful rerun flags:
  FORCE=1          rerun every step
  FORCE_KRAKEN=1   rerun Kraken2 only, while keeping earlier completed steps
EOF
}

require_command() {
  local command_name="$1"
  if ! command -v "$command_name" >/dev/null 2>&1; then
    echo "Missing command: $command_name" >&2
    exit 1
  fi
}

collect_kraken_targets() {
  local target_dir="$1"
  local -n targets_ref="$2"
  local candidate

  targets_ref=()

  shopt -s nullglob
  for candidate in "$target_dir"/combined_[0-9]*_unique.fasta; do
    targets_ref+=("$candidate")
  done
  shopt -u nullglob

  for candidate in \
    "$target_dir"/combined_consensussequences.fasta \
    "$target_dir"/combined_nogroup_unique.fasta; do
    if [[ -f "$candidate" ]]; then
      targets_ref+=("$candidate")
    fi
  done
}

kraken_combined_input_current() {
  local combined_input="$1"
  shift
  local file

  [[ -s "$combined_input" ]] || return 1

  for file in "$@"; do
    if [[ "$file" -nt "$combined_input" ]]; then
      return 1
    fi
  done

  return 0
}

kraken_outputs_current() {
  local kraken_output="$1"
  local kraken_report="$2"
  shift 2
  local file

  [[ -s "$kraken_output" && -s "$kraken_report" ]] || return 1

  for file in "$@"; do
    if [[ "$file" -nt "$kraken_output" || "$file" -nt "$kraken_report" ]]; then
      return 1
    fi
  done

  return 0
}

amplicon_sorter_outputs_exist() {
  local filtered_dir="$1"
  local target_dir="$filtered_dir/combined"
  local targets=()

  [[ -d "$target_dir" ]] || return 1

  collect_kraken_targets "$target_dir" targets
  if [[ "${#targets[@]}" -gt 0 ]]; then
    return 0
  fi

  [[ -s "$target_dir/results.txt" || -s "$filtered_dir/results.csv" ]]
}

validate_setup() {
  require_command fastplong
  require_command python
  require_command script

  if [[ "$KRAKEN_USE_DOCKER" == "1" ]]; then
    require_command docker
  else
    require_command kraken2
  fi

  [[ -f "$ADAPTER_FASTA" ]] || {
    echo "Missing adapter FASTA: $ADAPTER_FASTA" >&2
    exit 1
  }

  [[ -f "$AMPLICON_SORTER_SCRIPT" ]] || {
    echo "Missing amplicon sorter script: $AMPLICON_SORTER_SCRIPT" >&2
    exit 1
  }

  if [[ "$KRAKEN_USE_DOCKER" != "1" && ! -d "$KRAKEN_DB" ]]; then
    echo "Missing Kraken2 database directory: $KRAKEN_DB" >&2
    exit 1
  fi
}

filter_reads() {
  local folder="$1"
  local filtered_dir="$folder/filtered"
  local logs_dir="$filtered_dir/logs"
  local done_file="$filtered_dir/filtering.done"

  if [[ "$FORCE" != "1" && -f "$done_file" ]]; then
    echo "[$folder] filter_reads: already done"
    return
  fi

  mkdir -p "$logs_dir"

  shopt -s nullglob
  local files=("$folder"/*.fastq.gz)
  shopt -u nullglob

  if [[ "${#files[@]}" -eq 0 ]]; then
    echo "[$folder] no *.fastq.gz files found directly inside folder" >&2
    return 1
  fi

  local file filename base out html per_read_log
  for file in "${files[@]}"; do
    filename="$(basename "$file")"
    base="${filename%.fastq.gz}"
    out="$filtered_dir/${base}_filtered.fastq.gz"
    html="$filtered_dir/${base}_report.html"
    per_read_log="$logs_dir/${base}.log"

    if [[ "$FORCE" != "1" && -s "$out" ]]; then
      echo "[$folder] fastplong: skipping $filename"
      continue
    fi

# // TODO: check filtering options, find a way to print the barcodes trimmed sequences to log as currently is not appearing
    echo "[$folder] fastplong: $filename"

    local fastplong_cmd=(
      fastplong
      -i "$file"
      -o "$out"
      --unqualified_percent_limit "$FASTPLONG_UNQUALIFIED_PERCENT_LIMIT"
      --qualified_quality_phred "$FASTPLONG_QUALIFIED_QUALITY_PHRED"
      --mean_qual "$FASTPLONG_MEAN_QUAL"
      --length_required "$FASTPLONG_LENGTH_REQUIRED"
      --adapter_fasta "$ADAPTER_FASTA"
      -h "$html"
    )

    printf 'Command:' > "$per_read_log"
    printf ' %q' "${fastplong_cmd[@]}" >> "$per_read_log"
    printf '\n\n' >> "$per_read_log"

    # fastplong prints some adapter-trimming details only when attached to a TTY.
    # script provides a pseudo-terminal and records that screen output in the log.
    script -q -e -a -c "$(printf '%q ' "${fastplong_cmd[@]}")" "$per_read_log"
  done

  date > "$done_file"
}
# //TODO remove the .gz from combined.fastq.gz because the file is not compressed
combine_filtered_reads() {
  local folder="$1"
  local filtered_dir="$folder/filtered"
  local logs_dir="$filtered_dir/logs"
  local combined="$filtered_dir/combined.fastq.gz"
  local done_file="$filtered_dir/combined.done"

  if [[ "$FORCE" != "1" && -f "$done_file" && -s "$combined" ]]; then
    echo "[$folder] combine_filtered_reads: already done"
    return
  fi

  mkdir -p "$logs_dir"

  shopt -s nullglob
  local files=("$filtered_dir"/*_filtered.fastq.gz)
  shopt -u nullglob

  if [[ "${#files[@]}" -eq 0 ]]; then
    echo "[$folder] no filtered FASTQ files found in $filtered_dir" >&2
    return 1
  fi

  echo "[$folder] combining ${#files[@]} filtered FASTQ files"
  cat "${files[@]}" > "$combined"
  date > "$done_file"
}

run_amplicon_sorter() {
  local folder="$1"
  local filtered_dir="$folder/filtered"
  local combined="$filtered_dir/combined.fastq.gz"
  local logs_dir="$filtered_dir/logs"
  local log="$logs_dir/amplicon_sorter.log"
  local done_file="$filtered_dir/amplicon_sorter.done"

  if [[ "$FORCE" != "1" && -f "$done_file" ]]; then
    echo "[$folder] amplicon_sorter: already done"
    return
  fi

  if [[ "$FORCE" != "1" ]] && amplicon_sorter_outputs_exist "$filtered_dir"; then
    echo "[$folder] amplicon_sorter: outputs already exist"
    date > "$done_file"
    return
  fi

  [[ -s "$combined" ]] || {
    echo "[$folder] missing combined FASTQ: $combined" >&2
    return 1
  }

  mkdir -p "$logs_dir"

  echo "[$folder] amplicon_sorter"
  # Add -ra here if you want amplicon-sorter to reverse-complement reads.
  python "$AMPLICON_SORTER_SCRIPT" \
    -i "$combined" \
    -o "$filtered_dir/" \
    -min "$AMPLICON_MINLENGTH" \
    -max "$AMPLICON_MAXLENGTH" \
    -ar \
    -maxr "$AMPLICON_MAXREADS" \
    -np "$AMPLICON_THREADS" \
    -aln \
    -sg "$AMPLICON_SG" \
    -ss "$AMPLICON_SS" \
    -sc "$AMPLICON_SC" \
    $AMPLICON_EXTRA_ARGS \
    > "$log" 2>&1

  date > "$done_file"
}

run_kraken2() {
  local folder="$1"
  local target_dir="$folder/filtered/combined"
  local logs_dir="$folder/filtered/logs"
  local log="$logs_dir/kraken2.log"
  local done_file="$folder/filtered/kraken2.done"
  local combined_input="$target_dir/kraken2_combined_input.fasta"
  local kraken_output="$combined_input.kraken"
  local kraken_report="$combined_input.report"
  local memory_mapping=()

  [[ -d "$target_dir" ]] || {
    echo "[$folder] missing amplicon-sorter output directory: $target_dir" >&2
    return 1
  }

  mkdir -p "$logs_dir"
  : > "$log"

  if [[ "$KRAKEN_MEMORY_MAPPING" == "1" ]]; then
    memory_mapping=(--memory-mapping)
  fi

  local targets=()
  collect_kraken_targets "$target_dir" targets

  if [[ "${#targets[@]}" -eq 0 ]]; then
    echo "[$folder] no Kraken2 target FASTA files found in $target_dir" | tee -a "$log"
    date > "$done_file"
    return 0
  fi

  if [[ "$FORCE" != "1" && "$FORCE_KRAKEN" != "1" ]] \
    && kraken_outputs_current "$kraken_output" "$kraken_report" "${targets[@]}"; then
    echo "[$folder] kraken2: already done ($(basename "$kraken_output"))" | tee -a "$log"
    date > "$done_file"
    return 0
  fi

  if [[ "$FORCE" == "1" || "$FORCE_KRAKEN" == "1" ]] \
    || ! kraken_combined_input_current "$combined_input" "${targets[@]}"; then
    echo "[$folder] kraken2: combining ${#targets[@]} FASTA files into $(basename "$combined_input")" | tee -a "$log"
    cat "${targets[@]}" > "$combined_input"
  fi

  echo "[$folder] kraken2: $(basename "$combined_input")" | tee -a "$log"

  if [[ "$KRAKEN_USE_DOCKER" == "1" ]]; then
    local hostdir dbdir
    hostdir="$(realpath "$target_dir")"
    dbdir="$(realpath "$KRAKEN_DB")"
    docker run --rm \
      -v "$hostdir":/data \
      -v "$dbdir":/db \
      "$KRAKEN_DOCKER_IMAGE" kraken2 \
      --db /db \
      "${memory_mapping[@]}" \
      --confidence "$KRAKEN_CONFIDENCE" \
      --threads "$KRAKEN_THREADS" \
      "/data/$(basename "$combined_input")" \
      --output "/data/$(basename "$kraken_output")" \
      --report "/data/$(basename "$kraken_report")" \
      >> "$log" 2>&1
  else
    kraken2 \
      --db "$KRAKEN_DB" \
      "${memory_mapping[@]}" \
      --confidence "$KRAKEN_CONFIDENCE" \
      --threads "$KRAKEN_THREADS" \
      "$combined_input" \
      --output "$kraken_output" \
      --report "$kraken_report" \
      >> "$log" 2>&1
  fi

  date > "$done_file"
}

run_folder() {
  local folder="$1"

  [[ -d "$folder" ]] || {
    echo "Missing folder: $folder" >&2
    return 1
  }

  echo "[$folder] starting"
  filter_reads "$folder"
  combine_filtered_reads "$folder"
  run_amplicon_sorter "$folder"
  run_kraken2 "$folder"
  echo "[$folder] done"
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

  validate_setup

  local folder
  for folder in "$@"; do
    wait_for_slot
    run_folder "$folder" &
  done

  local failed=0
  for job in $(jobs -rp); do
    if ! wait "$job"; then
      failed=1
    fi
  done

  exit "$failed"
}

main "$@"
