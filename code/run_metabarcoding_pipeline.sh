#!/usr/bin/env bash
set -euo pipefail
CONFIG_FILE="${CONFIG_FILE:-config/pipeline.env}"
[[ ! -f "$CONFIG_FILE" ]] || source "$CONFIG_FILE"
: "${ADAPTER_FASTA:=data/ONT_24_native_barcode_sequences.fasta}"
: "${FASTPLONG_UNQUALIFIED_PERCENT_LIMIT:=50}"
: "${FASTPLONG_QUALIFIED_QUALITY_PHRED:=10}"
: "${FASTPLONG_MEAN_QUAL:=6}"
: "${FASTPLONG_LENGTH_REQUIRED:=300}"
: "${FASTPLONG_LENGTH_LIMIT:=1900}"
: "${FASTPLONG_TRIMMING_EXTENSION:=2}"
: "${SINTAX_REFERENCE:=data/sintax/MIDORI2_UNIQ_NUC_SP_GB272_CO1_SINTAX.fasta}"
: "${SINTAX_DB:=${SINTAX_REFERENCE%.fasta}.udb}"
: "${SINTAX_OUTPUT_ROOT:=data/sintax}"
: "${SINTAX_SUFFIX:=sintaxSP}"
: "${SINTAX_THREADS:=8}"
: "${SINTAX_CUTOFF:=1}"
: "${SINTAX_RANDSEED:=666}"
: "${SINTAX_MIN_LENGTH:=300}"
: "${SINTAX_MAX_LENGTH:=2000}"
: "${MAX_JOBS:=1}"
: "${FORCE:=0}"
: "${FORCE_SINTAX:=0}"

prepare_database() (
    # Called only when classification is needed, before starting folder workers.
    [[ ! -s "$SINTAX_DB" ]] || return 0
    command -v vsearch >/dev/null
    mkdir -p "$(dirname "$SINTAX_DB")" "$(dirname "$SINTAX_REFERENCE")"
    local temp_ref="" temp_db=""
    trap '[[ -z "$temp_ref" ]] || rm -f -- "$temp_ref"; [[ -z "$temp_db" ]] || rm -f -- "$temp_db"' EXIT
    if [[ ! -s "$SINTAX_REFERENCE" ]]; then
        [[ -s "${SINTAX_REFERENCE}.gz" ]] || { echo "Missing reference: $SINTAX_REFERENCE[.gz]" >&2; exit 1; }
        temp_ref="$(mktemp "${SINTAX_REFERENCE}.XXXXXX")"
        gzip -dc -- "${SINTAX_REFERENCE}.gz" > "$temp_ref"
        [[ -s "$temp_ref" ]]
        mv -- "$temp_ref" "$SINTAX_REFERENCE"
        temp_ref=""
    fi
    temp_db="$(mktemp "${SINTAX_DB}.XXXXXX")"
    vsearch --makeudb_usearch "$SINTAX_REFERENCE" --output "$temp_db"
    [[ -s "$temp_db" ]]
    mv -- "$temp_db" "$SINTAX_DB"
    temp_db=""
)

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
      --trimming_extension "$FASTPLONG_TRIMMING_EXTENSION"
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

run_folder() (
    local folder="$1" sample output tmp=""
    sample="$(basename "$folder")"
    output="$SINTAX_OUTPUT_ROOT/$sample/${sample}_${SINTAX_SUFFIX}.tsv"
    # Preserve completed results even if their parameters differ from today's settings.
    if [[ -f "$output" && "$FORCE" != 1 && "$FORCE_SINTAX" != 1 ]]; then
        echo "[$sample] existing SINTAX output; skipping: $output"
        return 0
    fi
    if [[ ! -s "$folder/filtered/asinput.fastq.gz" || "$FORCE" == 1 ]]; then
        command -v fastplong >/dev/null
        [[ -f "$ADAPTER_FASTA" ]] || { echo "Missing adapter FASTA: $ADAPTER_FASTA" >&2; exit 1; }
        filter_reads "$folder"
        combine_filtered_reads "$folder"
    fi
    mkdir -p "$(dirname "$output")"
    tmp="$(mktemp "${output}.XXXXXX")"
    trap '[[ -z "$tmp" ]] || rm -f -- "$tmp"' EXIT
    vsearch --sintax "$folder/filtered/asinput.fastq.gz" --db "$SINTAX_DB" \
        --sintax_random --strand both --sintax_cutoff "$SINTAX_CUTOFF" \
        --threads "$SINTAX_THREADS" --randseed "$SINTAX_RANDSEED" \
        --maxseqlength "$SINTAX_MAX_LENGTH" --minseqlength "$SINTAX_MIN_LENGTH" \
        --tabbedout "$tmp"
    mv -- "$tmp" "$output"
    tmp=""
)

main() {
    if [[ "$#" -eq 0 || "${1:-}" == --help ]]; then
        echo 'Usage: bash code/run_metabarcoding_pipeline.sh [--dry-run] FOLDER [FOLDER ...]'
        echo 'Activate taxonomy first. Each FOLDER is a sample folder; its basename identifies outputs.'
        return 0
    fi
    local dry_run=0
    if [[ "$1" == --dry-run ]]; then dry_run=1; shift; fi
    [[ "$#" -gt 0 && "$MAX_JOBS" =~ ^[1-9][0-9]*$ ]] || return 2
    local folder sample output needs_db=0 failed=0 pid
    local -a folders=() pids=()
    local -A seen=()
    for folder in "$@"; do
        [[ -d "$folder" ]] || { echo "Missing folder: $folder" >&2; return 1; }
        folder="$(realpath "$folder")"
        sample="$(basename "$folder")"
        [[ -z "${seen[$sample]:-}" ]] || { echo "Duplicate sample name: $sample; use separate output roots for separate datasets" >&2; return 1; }
        seen[$sample]=1
        folders+=("$folder")
        output="$SINTAX_OUTPUT_ROOT/$sample/${sample}_${SINTAX_SUFFIX}.tsv"
        if [[ -f "$output" && "$FORCE" != 1 && "$FORCE_SINTAX" != 1 ]]; then
            echo "[$sample] reuse $output"
        else
            needs_db=1
            echo "[$sample] filtered/asinput.fastq.gz -> $output (cutoff=$SINTAX_CUTOFF)"
        fi
    done
    [[ "$dry_run" != 1 ]] || return 0
    if [[ "$needs_db" == 1 ]]; then prepare_database; fi
    for folder in "${folders[@]}"; do
        run_folder "$folder" &
        pids+=("$!")
        if [[ "${#pids[@]}" -ge "$MAX_JOBS" ]]; then
            for pid in "${pids[@]}"; do wait "$pid" || failed=1; done
            pids=()
        fi
    done
    for pid in "${pids[@]}"; do wait "$pid" || failed=1; done
    return "$failed"
}
main "$@"
