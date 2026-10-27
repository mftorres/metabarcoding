# Export cleaned reads by best forward-primer tag

`code/export_primer_fasta.py` is a standalone exporter alongside
`code/export_sintax_fasta.py`. Selection starts from **best forward-tag hits on raw
reads**, while the output sequences come from **cleaned reads**. It does not rerun
primer detection, cleaning, SINTAX, notebooks or phylogenetics.

From the repository root:

```bash
conda activate taxonomy
python code/export_primer_fasta.py \
  --primer bc1087 \
  --barcode barcode15 \
  --analysis-dir data/sintax_analysis/mergedseq5_20261021 \
  --output data/primer_exports/barcode15_bc1087.fasta \
  --threads 4
```

Use a new output basename for every export. The script refuses to overwrite the
FASTA or any of its companion tables/logs.

## Selection and read tracing

1. Read the existing SeqKit-derived **best-hit table** recorded in the completed
   analysis. Match the dataset and requested barcode folder, and select only the
   `forward_patternName` field. Other/lower-ranked hits in the original locate
   output and reverse-primer hits do not select a read.
2. Map original names to canonical PCR tag names using the analysis's
   `primer_name_mapping.json`. For this dataset, `forward_tagged_primer_44` maps to
   **bc1087**. Already-canonical names are also accepted. If a saved mapping has
   several equivalent canonical names separated by semicolons, membership selects
   the read and the full mapping is retained in the audit table.
3. Use `seqkit grep --pattern-file` to locate the selected raw IDs in all top-level
   raw `.fastq`, `.fq`, `.fastq.gz` or `.fq.gz` files in the barcode folder. The raw
   filename is recorded. Missing raw IDs or duplicate raw IDs cause an error.
4. Find corresponding cleaned queries in `filtered/asinput.fastq.gz`, using the
   original read ID, the parent of `split-by-adapter-left-` / `split-by-adapter-right-`
   IDs, or an explicit `parent_read_id` in the FASTQ header. All surviving fragments
   of a selected parent are exported. A parent-level tag hit does **not** imply
   that every exported fragment itself contains the primer/tag sequence.
5. Annotate those queries with the **lowest trusted taxon** (lowest rank with
   confidence >=0.80) and its contamination flag from the completed analysis.
   Taxonomic confidence does not determine primer-based selection. Then extract the cleaned queries by
   exact ID with SeqKit and write annotated FASTA. Sequences remain in their
   cleaned orientation and are not trimmed again or reverse-complemented.

Reads removed by cleaning have no cleaned sequence to export; they are counted
and recorded as `no_clean_read`. A selected cleaned query without a SINTAX row
causes an error rather than receiving invented taxonomy. A valid but unassigned
SINTAX row is exported as a header with taxon `unassigned` and flag `nocontamination`, matching the
existing SINTAX exporter. Here `nocontamination` means the assigned label was not
flagged in the analysis, not proof that the read is uncontaminated.

## Headers and outputs

The FASTA header uses exactly the existing SINTAX export convention and sanitiser:

```text
>folder-forwardprimer-reverseprimer-lowest_trusted_taxon-contamination_or_nocontamination-cleaned_readID
```

For example:

```text
>barcode15-bc1087-bc1100rc-Tipula_maxima_123-contamination-split-by-adapter-left-READ_UUID
```

The header annotation is the **lowest trusted taxon**, the lowest rank with confidence **>=0.80**; a species at
0.79 falls back to a trusted genus/family/higher rank. No trusted rank gives
`unassigned`. Both detected primer names come from the selected raw parent's best
hit; absent forward/reverse hits are `none` (`unknown` in general taxon exports
when tracking is unavailable). `taxapass/taxafail` remains in audit tables but is
not included in the new six-field FASTA header. Spaces become underscores and
embedded hyphens in non-ID fields are escaped. Cleaned read IDs remain unchanged.

For `barcode15_bc1087.fasta`, the script also creates:

- `barcode15_bc1087.reads.tsv`: one row per exported query plus one row per selected
  raw parent with no surviving clean query. Columns: `source_folder`,
  `forward_primers_detected`, `reverse_primers_detected`, `raw_ont_read_id`, `raw_fastq_file`, `readid`,
  `clean_fastq_file`, `assigned_rank`, `assigned_label`, `assigned_confidence`,
  `taxonomy_status`, `contamination`, `export_header`, `status`.
- `barcode15_bc1087.counts.tsv`: per-folder and total counts, with columns
  `source_folder`, `raw_best_tag_hits`, `raw_parents_with_clean_reads`,
  `clean_reads_exported`, `raw_parents_without_clean_reads`. Split fragments can
  make the cleaned-read count exceed the number of raw parents.
- `barcode15_bc1087.log`: input provenance, SeqKit commands and extraction counts.

All tables are tab-delimited with headers. No best-tag hits produces an explicit
error. If there are raw hits but every one was removed by cleaning, the exporter
writes an empty FASTA with the explanatory counts/audit tables. Files are published
only after every requested folder validates.

## Optional unclassified folder

Add `--include-unclassified` to search both the requested barcode and the
`unclassified` folder. This is **off by default**. Folder selection remains separate:
`unclassified` sequences retain `unclassified` in their headers and audit rows.
The flag selects that folder's own best bc1087 hits; it does not relabel them as
barcode15 or search every other barcode.

Both folders must have raw files, best-hit tracking, cleaned reads and SINTAX
outputs included in the supplied completed analysis. The current
`mergedseq5_20261021` analysis covers barcode01–barcode21, so it cannot yet supply
annotated `unclassified` exports. The script reports this prerequisite without
launching classification or changing the analysis. After separately processing
unclassified and including it in a completed analysis, supply that analysis with
the same flag.

Optional overrides: `--primer-best-hits PATH` selects another compatible best-hit
TSV and `--merged-dir PATH` changes the raw/clean dataset root. Normally the saved
analysis paths should be used. The script checks the saved SINTAX size/mtime
against `input_manifest.tsv` and rejects changed classifier outputs.

## Tests

```bash
conda activate taxonomy
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s code -p test_export_primer_fasta.py -v
```

These use real SeqKit with temporary raw/clean/SINTAX fixtures. They verify parent
ID and adapter-split tracing, forward-only canonical tag selection, identical taxon
export headers, contamination and unassigned labels, optional unclassified output,
missing-raw-ID failures, and refusal to overwrite outputs. Existing data are untouched.

## Consistent selection and annotation terminology

This exporter selects reads by their raw parent's **best forward-primer hit**.
The taxon exporter instead selects reads by a **reported taxon label** in the full
SINTAX lineage with confidence >=0.80 on the queried label itself. Both write the **lowest
trusted taxon** (confidence >=0.80) in the sequence header, or `unassigned` if no
rank qualifies. The header taxon can therefore differ from the queried taxon.
Neither exporter changes the confidence-based species-summary or contamination rules.
