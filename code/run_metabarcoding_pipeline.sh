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
: "${FASTPLONG_LENGTH_REQUIRED:=300}"
: "${FASTPLONG_LENGTH_LIMIT:=1500}"
: "${AMPLICON_SORTER_SCRIPT:=code/amplicon_sorter_mod.py}"
: "${AMPLICON_READ_MAP_SCRIPT:=code/build_amplicon_kraken_read_map.py}"
: "${AMPLICON_THREADS:=6}"
: "${AMPLICON_MINLENGTH:=300}"
: "${AMPLICON_MAXLENGTH:=100000000}"
: "${AMPLICON_MAXREADS:=10000000000000000}"
: "${AMPLICON_SG:=99.0}"
: "${AMPLICON_SS:=99.0}"
: "${AMPLICON_SC:=99.0}"
: "${AMPLICON_EXTRA_ARGS:=}"
: "${AMPLICON_TMP_ROOT:=/tmp/metabarcoding_amplicon_sorter}"
: "${NONUMT_MIN_AA_LENGTH:=200}"
: "${NONUMT_TRANSLATION_TABLE:=5}"
: "${NONUMT_SEQKIT_THREADS:=6}"
: "${SEQKIT_CONDA_ENV:=umitracker}"
: "${KRAKEN_THREADS:=4}"
: "${KRAKEN_DB:=/mnt/d/kraken2DB/k2_core_nt_2025}"
: "${KRAKEN_CONFIDENCE:=0}"
: "${KRAKEN_MEMORY_MAPPING:=1}"
: "${KRAKEN_USE_DOCKER:=0}"
: "${KRAKEN_DOCKER_IMAGE:=staphb/kraken2}"
: "${MAX_JOBS:=1}"
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

set_seqkit_command() {
  if command -v seqkit >/dev/null 2>&1; then
    SEQKIT_CMD=(seqkit)
    return
  fi

  if command -v conda >/dev/null 2>&1; then
    SEQKIT_CMD=(conda run -n "$SEQKIT_CONDA_ENV" seqkit)
    return
  fi

  echo "Missing seqkit. Activate the ${SEQKIT_CONDA_ENV} conda environment or make conda available." >&2
  exit 1
}

collect_kraken_target() {
  local target_dir="$1"
  local -n target_ref="$2"
  local candidate
  local targets=()

  target_ref=""

  shopt -s nullglob
  for candidate in "$target_dir"/*consensussequences_noNUMT_"${NONUMT_MIN_AA_LENGTH}"filter.fasta; do
    [[ -s "$candidate" ]] || continue
    targets+=("$candidate")
  done
  shopt -u nullglob

  if [[ "${#targets[@]}" -gt 1 ]]; then
    echo "Expected one noNUMT consensus FASTA in $target_dir, found ${#targets[@]}" >&2
    return 1
  fi

  if [[ "${#targets[@]}" -eq 1 ]]; then
    target_ref="${targets[0]}"
  fi
}

kraken_read_map_current() {
  local read_map="$1"
  shift
  local file

  [[ -s "$read_map" ]] || return 1

  for file in "$@"; do
    if [[ "$file" -nt "$read_map" ]]; then
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

  [[ -d "$target_dir" ]] || return 1

  if [[ -s "$target_dir/asinput_consensussequences.fasta" ]]; then
    [[ -s "$target_dir/amplicon_sorter_read_map.tsv" ]] || return 1
    return 0
  fi

  [[ -s "$target_dir/results.txt" || -s "$filtered_dir/results.csv" ]]
}

amplicon_sorter_output_usable() {
  local target_dir="$1"
  local log="$2"

  [[ -d "$target_dir" && -s "$target_dir/results.txt" ]] || return 1

  if [[ -s "$target_dir/asinput_consensussequences.fasta" ]]; then
    [[ -s "$target_dir/amplicon_sorter_read_map.tsv" ]] || return 1
    return 0
  fi

  # amplicon_sorter exits cleanly for samples with too few usable reads.
  # Treat that as a completed no-target run, but reject other no-FASTA outputs.
  [[ -f "$log" ]] && grep -q "Number of usable sequences is lower than 5, quitting job" "$log"
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
  [[ -f "$AMPLICON_READ_MAP_SCRIPT" ]] || {
    echo "Missing amplicon/Kraken read-map script: $AMPLICON_READ_MAP_SCRIPT" >&2
    exit 1
  }
  set_seqkit_command

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
      --length_limit "$FASTPLONG_LENGTH_LIMIT"
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
combine_filtered_reads() {
  local folder="$1"
  local filtered_dir="$folder/filtered"
  local logs_dir="$filtered_dir/logs"
  local combined="$filtered_dir/asinput.fastq.gz"
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
  local combined="$filtered_dir/asinput.fastq.gz"
  local target_dir="$filtered_dir/combined"
  local logs_dir="$filtered_dir/logs"
  local log="$logs_dir/amplicon_sorter.log"
  local stdout_log="$logs_dir/amplicon_sorter.stdout.log"
  local stderr_log="$logs_dir/amplicon_sorter.stderr.log"
  local done_file="$filtered_dir/amplicon_sorter.done"
  local sample_name tmp_parent tmp_target staging root_file

  mkdir -p "$logs_dir"

  if [[ "$FORCE" != "1" && -f "$done_file" ]] && amplicon_sorter_output_usable "$target_dir" "$log"; then
    echo "[$folder] amplicon_sorter: already done"
    return
  fi

  if [[ "$FORCE" != "1" ]] && amplicon_sorter_output_usable "$target_dir" "$log"; then
    echo "[$folder] amplicon_sorter: outputs already exist"
    date > "$done_file"
    return
  fi

  [[ -s "$combined" ]] || {
    echo "[$folder] missing combined FASTQ: $combined" >&2
    return 1
  }

  mkdir -p "$AMPLICON_TMP_ROOT"
  sample_name="$(basename "$folder")"
  tmp_parent="$(mktemp -d "$AMPLICON_TMP_ROOT/${sample_name}.amplicon_sorter.XXXXXX")"
  tmp_target="$tmp_parent/asinput"
  staging="$filtered_dir/.combined.copying.$$"

  echo "[$folder] amplicon_sorter: writing temporary output to $tmp_parent"
  # Add -ra here if you want amplicon-sorter to reverse-complement reads.
  # Running amplicon_sorter using a combined file with all filtered reads from the same barcode
  # It is less efficient for storage but amplicon_sorter will grab any fastq, fastq, gz file even if not relevant for the analyses
  # need to change the folder structure in the future, to avoid additional read files and be able to pass a folder as input
  # or maybe keeping all filtered reads in a single file and pass it as input for amplicon_sorter is better
  # so that all reads can be compared to all
  if ! python "$AMPLICON_SORTER_SCRIPT" \
    -i "$combined" \
    -o "$tmp_parent/" \
    -min "$AMPLICON_MINLENGTH" \
    -max "$AMPLICON_MAXLENGTH" \
    -ar \
    -maxr "$AMPLICON_MAXREADS" \
    -np "$AMPLICON_THREADS" \
    -sg "$AMPLICON_SG" \
    -ss "$AMPLICON_SS" \
    -sc "$AMPLICON_SC" \
    $AMPLICON_EXTRA_ARGS \
    > "$stdout_log" 2> "$stderr_log"; then
    cat "$stdout_log" "$stderr_log" > "$log"
    echo "[$folder] amplicon_sorter failed; stdout: $stdout_log; stderr: $stderr_log" | tee -a "$log" >&2
    echo "[$folder] keeping temporary output for inspection: $tmp_parent" | tee -a "$log" >&2
    return 1
  fi

  cat "$stdout_log" "$stderr_log" > "$log"

  if ! amplicon_sorter_output_usable "$tmp_target" "$log"; then
    echo "[$folder] amplicon_sorter did not produce usable output in $tmp_target" | tee -a "$log" >&2
    echo "[$folder] keeping temporary output for inspection: $tmp_parent" | tee -a "$log" >&2
    return 1
  fi

  rm -rf "$staging"
  cp -a "$tmp_target" "$staging"

  if ! amplicon_sorter_output_usable "$staging" "$log"; then
    echo "[$folder] copied amplicon_sorter output is not usable: $staging" | tee -a "$log" >&2
    rm -rf "$staging"
    return 1
  fi

  rm -rf "$target_dir"
  mv "$staging" "$target_dir"

  shopt -s nullglob
  for root_file in "$tmp_parent"/*; do
    if [[ -f "$root_file" ]]; then
      cp -a "$root_file" "$filtered_dir/"
    fi
  done
  shopt -u nullglob

  if ! amplicon_sorter_output_usable "$target_dir" "$log"; then
    echo "[$folder] final amplicon_sorter output is not usable after copy: $target_dir" | tee -a "$log" >&2
    return 1
  fi

  rm -rf "$tmp_parent"
  echo "[$folder] amplicon_sorter: copied output to $target_dir and removed temporary output"
  date > "$done_file"
}

run_kraken2() {
  local folder="$1"
  local filtered_dir="$folder/filtered"
  local target_dir="$folder/filtered/combined"
  local logs_dir="$folder/filtered/logs"
  local log="$logs_dir/kraken2.log"
  local done_file="$folder/filtered/kraken2.done"
  local kraken_input=""
  local read_map
  local kraken_output
  local kraken_report
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

  collect_kraken_target "$target_dir" kraken_input

  if [[ -z "$kraken_input" ]]; then
    echo "[$folder] no Kraken2 target FASTA files found in $target_dir" | tee -a "$log"
    rm -f "$done_file"
    return 1
  fi

  read_map="${kraken_input%.fasta}.read_map.tsv"
  kraken_output="$kraken_input.kraken"
  kraken_report="$kraken_input.report"

  if [[ "$FORCE" != "1" && "$FORCE_KRAKEN" != "1" ]] \
    && kraken_outputs_current "$kraken_output" "$kraken_report" "$kraken_input"; then
    echo "[$folder] kraken2: already done ($(basename "$kraken_output"))" | tee -a "$log"
    date > "$done_file"
    return 0
  fi

  if [[ "$FORCE" == "1" || "$FORCE_KRAKEN" == "1" ]] \
    || ! kraken_read_map_current "$read_map" "$kraken_input"; then
    echo "[$folder] kraken2: writing read map $(basename "$read_map")" | tee -a "$log"
    python "$AMPLICON_READ_MAP_SCRIPT" \
      --output "$read_map" \
      --targets "$kraken_input" \
      >> "$log" 2>&1
  fi

  echo "[$folder] kraken2: $(basename "$kraken_input")" | tee -a "$log"

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
      "/data/$(basename "$kraken_input")" \
      --output "/data/$(basename "$kraken_output")" \
      --report "/data/$(basename "$kraken_report")" \
      >> "$log" 2>&1
  else
    kraken2 \
      --db "$KRAKEN_DB" \
      "${memory_mapping[@]}" \
      --confidence "$KRAKEN_CONFIDENCE" \
      --threads "$KRAKEN_THREADS" \
      "$kraken_input" \
      --output "$kraken_output" \
      --report "$kraken_report" \
      >> "$log" 2>&1
  fi

  date > "$done_file"
}

run_consensus_noNUMT_filter() {
  local folder="$1"
  local combined_dir="$folder/filtered/combined"
  local logs_dir="$folder/filtered/logs"
  local log="$logs_dir/noNUMT_filter.log"
  local consensus_fasta
  local input_name
  local input_stem
  local aa_output_fasta
  local nt_output_fasta
  local passing_headers
  local input_sequence_count
  local aa_output_sequence_count
  local nt_output_sequence_count
  local output_count=0

  mkdir -p "$logs_dir"
  echo "[$folder] noNUMT: filtering consensus sequences" | tee "$log"

  [[ -d "$combined_dir" ]] || {
    echo "[$folder] noNUMT: missing combined folder: $combined_dir" | tee -a "$log" >&2
    return 1
  }

  shopt -s nullglob
  local consensus_fastas=("$combined_dir"/*consensussequences.fasta)
  shopt -u nullglob

  if [[ "${#consensus_fastas[@]}" -eq 0 ]]; then
    echo "[$folder] noNUMT: no *consensussequences.fasta files found" | tee -a "$log"
    return 0
  fi

  for consensus_fasta in "${consensus_fastas[@]}"; do
    [[ "$consensus_fasta" != *noNUMT* ]] || continue

    input_name="$(basename "$consensus_fasta")"
    input_stem="${input_name%.fasta}"
    aa_output_fasta="$combined_dir/${input_stem}_aa_noNUMT_${NONUMT_MIN_AA_LENGTH}filter.fasta"
    nt_output_fasta="$combined_dir/${input_stem}_noNUMT_${NONUMT_MIN_AA_LENGTH}filter.fasta"
    input_sequence_count="$(grep -c '^>' "$consensus_fasta" || true)"

    if [[ "$FORCE" != "1" && -s "$nt_output_fasta" && "$nt_output_fasta" -nt "$consensus_fasta" ]]; then
      echo "[$folder] noNUMT: nucleotide output already exists: $nt_output_fasta" | tee -a "$log"
      aa_output_sequence_count="$(grep -c '^>' "$aa_output_fasta" 2>/dev/null || true)"
      nt_output_sequence_count="$(grep -c '^>' "$nt_output_fasta" || true)"
      echo "[$folder] noNUMT counts: input_nt=$input_sequence_count passing_aa_frames=$aa_output_sequence_count output_nt_noNUMT=$nt_output_sequence_count" | tee -a "$log"
      output_count=$((output_count + 1))
      continue
    fi

    echo "[$folder] noNUMT: translating and filtering $consensus_fasta" | tee -a "$log"
    echo "[$folder] noNUMT: writing amino-acid diagnostic output $aa_output_fasta" | tee -a "$log"
    echo "[$folder] noNUMT: writing cleaned nucleotide Kraken input $nt_output_fasta" | tee -a "$log"
    passing_headers="$(mktemp)"

    # Translate in all six frames and keep only amino-acid records at or above the
    # noNUMT length threshold. This amino-acid FASTA is diagnostic/filtering
    # output only; Kraken2 must receive nucleotide sequences.
    "${SEQKIT_CMD[@]}" translate "$consensus_fasta" \
      --frame 6 \
      --append-frame \
      -T "$NONUMT_TRANSLATION_TABLE" \
      -j "$NONUMT_SEQKIT_THREADS" \
      --trim \
      2>> "$log" \
      | awk -v min_len="$NONUMT_MIN_AA_LENGTH" '
        function flush_record() {
          if (header != "" && length(sequence) >= min_len) {
            print header
            for (i = 1; i <= length(sequence); i += 60) {
              print substr(sequence, i, 60)}
          }}
        /^>/ {flush_record(); header = $0; sequence = ""; next}
        {sequence = sequence $0}
        END {flush_record()}
      ' > "$aa_output_fasta"

    # Strip the frame suffix from passing amino-acid headers, deduplicate the
    # original consensus IDs, and fetch those IDs from the nucleotide consensus
    # FASTA. The resulting nucleotide FASTA is the Kraken2 input.
    awk '
      /^>/ {
        header = $0
        sub(/^>/, "", header)
        sub(/\r$/, "", header)
        sub(/[[:space:]]+$/, "", header)
        sub(/[_[:space:]]frame=-?[0-9]+.*/, "", header)
        if (!(header in seen)) {
          seen[header] = 1
          print header
        }
      }
    ' "$aa_output_fasta" > "$passing_headers"

    awk '
      FILENAME == ARGV[1] {
        keep[$0] = 1
        next
      }
      /^>/ {
        header = $0
        sub(/^>/, "", header)
        sub(/\r$/, "", header)
        sub(/[[:space:]]+$/, "", header)
        sub(/[_[:space:]]frame=-?[0-9]+.*/, "", header)
        write_record = (header in keep)
        if (write_record) {
          print ">" header
        }
        next
      }
      write_record {
        print
      }
    ' "$passing_headers" "$consensus_fasta" > "$nt_output_fasta"

    rm "$passing_headers"

    aa_output_sequence_count="$(grep -c '^>' "$aa_output_fasta" || true)"
    nt_output_sequence_count="$(grep -c '^>' "$nt_output_fasta" || true)"
    echo "[$folder] noNUMT counts: input_nt=$input_sequence_count passing_aa_frames=$aa_output_sequence_count output_nt_noNUMT=$nt_output_sequence_count" | tee -a "$log"

    output_count=$((output_count + 1))
  done

  echo "[$folder] noNUMT: nucleotide FASTA files written or already present: $output_count" | tee -a "$log"
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
  run_consensus_noNUMT_filter "$folder"
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
