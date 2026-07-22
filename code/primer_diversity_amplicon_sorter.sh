#!/usr/bin/env bash
set -euo pipefail

MERGED_DIR="${1:-data/mergedseq5}"
MERGE_NAME="$(basename "$MERGED_DIR")"
PRIMER_BEST_HITS="${2:-}"
OUTPUT_FILE="${3:-${MERGED_DIR}/${MERGE_NAME}_amplicon_sorter_primer_diversity.tsv}"

if [[ -z "$PRIMER_BEST_HITS" ]]; then
    PRIMER_BEST_HITS="$(find data/primertrack \
        -maxdepth 1 \
        -type f \
        -name "${MERGE_NAME}_seqkit_primer_umi_best_hits_*.txt" \
        | sort \
        | tail -n 1)"
fi

if [[ -z "$PRIMER_BEST_HITS" || ! -f "$PRIMER_BEST_HITS" ]]; then
    echo "Could not find primer best-hit table for ${MERGE_NAME}" >&2
    echo "Provide it as argument 2, for example:" >&2
    echo "  $0 ${MERGED_DIR} data/primertrack/${MERGE_NAME}_seqkit_primer_umi_best_hits_YYYYMMDD.txt" >&2
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT_FILE")"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT
PRIMER_BY_BARCODE_DIR="$TMP_DIR/primer_by_barcode"
mkdir -p "$PRIMER_BY_BARCODE_DIR"

printf 'barcode\tasinput_file\tcluster_group\tcluster_number\tgroup_number\tsequence_count\traw_read_count\tforward_detected_read_count\treverse_detected_read_count\tunique_forward_primer_count\tunique_reverse_primer_count\tforward_primers_detected\treverse_primers_detected\tforward_primer_counts\treverse_primer_counts\tfile_path\n' > "$OUTPUT_FILE"

awk -F '\t' -v merge_name="$MERGE_NAME" -v primer_dir="$PRIMER_BY_BARCODE_DIR" '
    FNR == 1 {
        for (i = 1; i <= NF; i++) {
            primer_col[$i] = i
        }
        next
    }
    $primer_col["merged_folder"] == merge_name {
        seq_id = $primer_col["seqID"]
        barcode = $primer_col["barcode_folder"]
        forward = $primer_col["forward_patternName"]
        reverse = $primer_col["reverse_patternName"]
        print seq_id "\t" forward "\t" reverse >> primer_dir "/" barcode "_primer_hits.tsv"
    }
' "$PRIMER_BEST_HITS"

for barcode_dir in "$MERGED_DIR"/barcode*; do
    [[ -d "$barcode_dir/filtered/combined" ]] || continue

    barcode="$(basename "$barcode_dir")"
    combined_dir="$barcode_dir/filtered/combined"
    amplicon_map="$combined_dir/amplicon_sorter_read_map.tsv"
    membership_file="$TMP_DIR/${barcode}_asinput_membership.tsv"
    barcode_primer_hits="$PRIMER_BY_BARCODE_DIR/${barcode}_primer_hits.tsv"

    if [[ ! -f "$amplicon_map" ]]; then
        echo "Skipping ${barcode}: missing ${amplicon_map}" >&2
        continue
    fi
    if [[ ! -f "$barcode_primer_hits" ]]; then
        : > "$barcode_primer_hits"
    fi

    : > "$membership_file"
    find "$combined_dir" \
        -maxdepth 1 \
        -type f \
        -name 'asinput_*_*.fasta' \
        -print0 \
        | sort -z \
        | while IFS= read -r -d '' fasta_file; do
            asinput_file="$(basename "$fasta_file")"
            awk -v asinput_file="$asinput_file" -v fasta_file="$fasta_file" '
                /^>/ {
                    read_id = $0
                    sub(/^>/, "", read_id)
                    sub(/[[:space:]].*/, "", read_id)
                    print asinput_file "\t" read_id "\t" fasta_file
                }
            ' "$fasta_file" >> "$membership_file"
        done

    awk -F '\t' -v barcode="$barcode" '
        FILENAME == ARGV[1] {
            seq_id = $1
            forward_primer[seq_id] = $2
            reverse_primer[seq_id] = $3
            next
        }
        FILENAME == ARGV[2] {
            if (FNR == 1) {
                for (i = 1; i <= NF; i++) {
                    map_col[$i] = i
                }
                next
            }
            internal_id = $map_col["amplicon_sorter_read_number"]
            ont_read_id = $map_col["ont_read_id"]
            internal_to_ont[internal_id] = ont_read_id
            next
        }
        FILENAME == ARGV[3] {
            asinput_file = $1
            internal_id = $2
            file_path = $3
            asinput_files[asinput_file] = 1
            file_paths[asinput_file] = file_path
            sequence_count[asinput_file] += 1

            ont_read_id = internal_to_ont[internal_id]
            if (ont_read_id == "") {
                next
            }

            if (raw_seen[asinput_file, ont_read_id] != 1) {
                raw_seen[asinput_file, ont_read_id] = 1
                raw_read_count[asinput_file] += 1
            }

            forward = forward_primer[ont_read_id]
            reverse = reverse_primer[ont_read_id]

            if (forward != "" && forward_seen_by_read[asinput_file, ont_read_id] != 1) {
                forward_seen_by_read[asinput_file, ont_read_id] = 1
                forward_detected_read_count[asinput_file] += 1
            }
            if (reverse != "" && reverse_seen_by_read[asinput_file, ont_read_id] != 1) {
                reverse_seen_by_read[asinput_file, ont_read_id] = 1
                reverse_detected_read_count[asinput_file] += 1
            }
            if (forward != "") {
                forward_count[asinput_file, forward] += 1
                if (forward_seen[asinput_file, forward] != 1) {
                    forward_seen[asinput_file, forward] = 1
                    forward_order[asinput_file] = forward_order[asinput_file] forward ";"
                    unique_forward_primer_count[asinput_file] += 1
                }
            }
            if (reverse != "") {
                reverse_count[asinput_file, reverse] += 1
                if (reverse_seen[asinput_file, reverse] != 1) {
                    reverse_seen[asinput_file, reverse] = 1
                    reverse_order[asinput_file] = reverse_order[asinput_file] reverse ";"
                    unique_reverse_primer_count[asinput_file] += 1
                }
            }
        }
        END {
            for (asinput_file in asinput_files) {
                stem = asinput_file
                sub(/\.fasta$/, "", stem)
                cluster_group = stem
                cluster_number = "none"
                group_number = "none"
                if (stem ~ /^asinput_.*_unique$/) {
                    cluster_number = stem
                    sub(/^asinput_/, "", cluster_number)
                    sub(/_unique$/, "", cluster_number)
                    group_number = "unique"
                    cluster_group = cluster_number "_unique"
                } else if (stem ~ /^asinput_[^_]+_[^_]+$/) {
                    cluster_group = stem
                    sub(/^asinput_/, "", cluster_group)
                    split(cluster_group, parts, "_")
                    cluster_number = parts[1]
                    group_number = parts[2]
                }

                forward_list = forward_order[asinput_file]
                reverse_list = reverse_order[asinput_file]
                sub(/;$/, "", forward_list)
                sub(/;$/, "", reverse_list)
                if (forward_list == "") {
                    forward_list = "none"
                }
                if (reverse_list == "") {
                    reverse_list = "none"
                }

                forward_counts = ""
                split(forward_order[asinput_file], forward_values, ";")
                for (forward_index in forward_values) {
                    forward = forward_values[forward_index]
                    if (forward != "") {
                        forward_counts = forward_counts forward ":" forward_count[asinput_file, forward] ";"
                    }
                }
                sub(/;$/, "", forward_counts)
                if (forward_counts == "") {
                    forward_counts = "none"
                }

                reverse_counts = ""
                split(reverse_order[asinput_file], reverse_values, ";")
                for (reverse_index in reverse_values) {
                    reverse = reverse_values[reverse_index]
                    if (reverse != "") {
                        reverse_counts = reverse_counts reverse ":" reverse_count[asinput_file, reverse] ";"
                    }
                }
                sub(/;$/, "", reverse_counts)
                if (reverse_counts == "") {
                    reverse_counts = "none"
                }

                printf "%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\n",
                    barcode,
                    asinput_file,
                    cluster_group,
                    cluster_number,
                    group_number,
                    sequence_count[asinput_file] + 0,
                    raw_read_count[asinput_file] + 0,
                    forward_detected_read_count[asinput_file] + 0,
                    reverse_detected_read_count[asinput_file] + 0,
                    unique_forward_primer_count[asinput_file] + 0,
                    unique_reverse_primer_count[asinput_file] + 0,
                    forward_list,
                    reverse_list,
                    forward_counts,
                    reverse_counts,
                    file_paths[asinput_file]
            }
        }
    ' "$barcode_primer_hits" "$amplicon_map" "$membership_file" \
        | sort -t $'\t' -k1,1 -k2,2 >> "$OUTPUT_FILE"
done

echo "Wrote $OUTPUT_FILE"
