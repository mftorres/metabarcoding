# Module 1: Cleaning and taxonomic assignment

The original pipeline is the cleaning and taxonomic-assignment module:
`code/run_metabarcoding_pipeline.sh` cleans reads and runs SINTAX, and the notebook
runner below generates the read-level analyses. The separate
[Module 2: phylogenetics](post_processing.md) uses its cleaned reads and completed
analysis to retrieve references, cluster with VSEARCH at 98% (local barcode/primer/taxon
groups separately; references globally), filter translations, align and build trees.
Module 2 also generates Baltic midpoint-rooted divergence-tree figures with
branch supports, organism/accession or lowest-trusted-taxon/readID labels, and a scale bar.
Mitogenomes are excluded from Module 2 references by default; its
`--include-mitogenomes` option allows them.

Activate `conda activate taxonomy` and run commands from the project root.
The current analysis uses cleaned reads at
`data/mergedseq5/<folder>/filtered/asinput.fastq.gz`, existing VSEARCH SINTAX results,
and existing primer-tracking best hits. No consensus generation, nucleotide-to-amino-acid
translation, or noNUMT filtering is performed in Module 1. Those translation and
filtering steps belong to the optional phylogenetics module.

## Run the post-classification notebooks

```bash
conda activate taxonomy
python code/run_sintax_notebooks.py \
  --output-dir data/sintax_analysis/my_new_run \
  --primer-hits data/primertrack/mergedseq5_seqkit_primer_umi_best_hits_20260731.txt
```

By default the requested folders are barcode01–barcode24. Use `--folders barcode01
barcode02` for a subset. Missing classifier outputs are recorded in
`input_manifest.tsv` and excluded from the analysis and PCR-plate denominator.
The runner never launches SINTAX, filtering or primer detection. A new output
folder is required. `--resume` resumes reporting only after parsing has completed.
`bash code/run_full_downstream_pipeline.sh` forwards the same command-line options.

The runner uses the Python environment from which it is launched for the notebook
kernel. It executes:

1. `code/sintax_output_parser.ipynb`: parses the available SINTAX files, joins primer
   results, calculates cross-plate prevalence, and writes all read/species tables.
2. `code/evaluate_tagtracking_taxid_reduced.ipynb` for each analysed folder: reviews
   read IDs, species and primer evidence and plots confident species read counts.
3. `code/sankey_target_species_detection.ipynb`: the retained filename now contains
   the across-folder read/species and contamination summaries, rather than the old
   consensus-derived target-species analysis.

Executed notebooks are saved under `code/executed_notebooks/sintax_<timestamp>/`.
`analysis.complete` marks successful table generation; `notebooks.complete` records
where the successful executed notebooks were saved. A failed or partial analysis
must not be used for FASTA exports.

For manual execution set `SINTAX_ANALYSIS_ROOT` to the new output directory before
opening the parser notebook. The per-folder notebook uses `SELECTED_BARCODE`.
The original `data/sintax/sintax_output_parser.ipynb` is also updated. Run only one
copy of the parser per output directory. The Tipuloidea notebook reads the combined
species summary and filters its four family names without creating new assignments.

## Confidence and taxonomy

Each read table retains the full original query header, the cleaned query `readid`,
original `raw_ont_read_id`, strand, full candidate taxonomy, and original column-four
`taxlabel`. Every rank (domain through strain where provided) has its candidate
label and `<rank>_bootstrap` confidence. Accepted rank labels and their confidences
are reported in separate `accepted_<rank>` and `accepted_<rank>_confidence` columns.

The analyses recompute acceptance at **confidence >=0.8** from the recorded values
in parentheses. The deepest qualifying rank supplies `assigned_rank_code`,
`assigned_rank`, `assigned_label`, `assigned_name`, `assigned_taxid`, and
`assigned_confidence`. Thus an old cutoff-1 column-four label can differ from the
new analysis assignment. The original label is preserved. These are rounded
bootstrap confidences reported by SINTAX, not calibrated probabilities.

`taxonomy_status` is `taxapass` only when the lowest accepted rank is genus or
species, otherwise `taxafail`. Unassigned reads remain in read tables with blank
assignment fields. Raw-parent metadata (or adapter-split prefixes) links fragments
to raw-read primer results; sequences are never taken from raw FASTQ for FASTA export.

## PCR-plate prevalence / contamination screen

The plate mapping is `saproxDNA/species_info_seq5/seq5_pcr_to_pools.tsv`.
Each taxon at each rank is present on a plate when at least one cleaned query
has confidence >=0.8 for that taxon/rank. Several barcodes mapping to one plate
count as one plate. The denominator is the number of distinct `PCR_plate` codes
among **analysed** folders, including any control codes if those folders are analysed.
Missing SINTAX files do not count as absence. Zero-call analysed samples do remain
in the denominator.

A taxon/rank receives `contamination` when it appears in at least 80% of those
plate codes; otherwise it receives `nocontamination`. For 21 plates, at least 17
must contain it. `taxon_plate_prevalence.tsv` records the count, denominator,
fraction, threshold, and supporting plate codes for every confidently observed taxon.
This is a prevalence screen, not proof of contamination. Broad ranks can legitimately
occur on almost every plate.

The read's main `contamination` flag applies to its **lowest accepted label**.
Separate `<rank>_contamination` columns preserve flags for all accepted ranks.
The species-summary flag applies to the species itself: an insect class being
widespread does not automatically mark every insect species as contamination.

## Tab-delimited outputs

Every table is a tab-delimited `.tsv` with column names in its first row. Exports explicitly enable headers. `TABLE_HEADERS.md` lists the ordered columns for every output type, and the notebook runner validates each header. JSON is used for run configuration and the primer-name mapping dictionary.
Per folder:

- `<folder>/<folder>_reads.tsv`: one row per cleaned SINTAX query, all rank/confidence
  fields, lowest accepted assignment, plate, mapped forward/reverse PCR primer names,
  source folder/file, original parent read ID, and contamination flags.
- `<folder>/<folder>_species_summary.tsv`: one row per unique candidate species
  label, with species name/taxid, genus, family, order, class, contamination flag,
  `median_confidence`, `read_count`, `raw_read_count`, `forward_primers_detected`, and `reverse_primers_detected`.
- `<folder>/<folder>_primer_taxonomy.tsv`: read-level accepted taxonomy grouped by
  detected forward/reverse PCR primer pair, with folder and plate provenance;
  reverse-only detections are retained.

`read_count` counts unique **cleaned query IDs** having species confidence >=0.8;
adapter-split fragments are separate cleaned queries. `raw_read_count` additionally
reports distinct original parent IDs. `median_confidence` is the median species
confidence among those passing queries. Candidate species without any passing
query remain visible with count zero and a blank confident median.
`candidate_read_count` and `median_candidate_confidence` describe all candidate calls.
Genus/family/order/class come from those candidate lineages; inconsistent labels for a
species are retained as semicolon-separated lists rather than silently discarded.
Forward/reverse-primer lists use the passing species reads, match raw parent IDs within
the dataset/barcode, and translate MFT names to PCR names by exact primer-sequence
identity. Only the mapped dictionary values are exported in biological tables.
The original-name columns (`forward_primer_names_original`, `forward_patternName`,
and `reverse_patternName`) are not included in these outputs. The input primer-tracking
tables are preserved, and `primer_name_mapping.json` retains the dictionary for audit.
For example, `forward_tagged_primer_1` maps to `bc1002`, and
`reverse_tagged_primer_1` maps to `bc1100rc`. The `bc...` names are the PCR names
(dictionary values), not the original MFT names (dictionary keys).
Read/primer-pair tables use `forward_PCR_primer` and `reverse_PCR_primer`;
aggregated tables use `forward_primers_detected` and `reverse_primers_detected`.
Lists are unique, semicolon-separated names; unknown names cause an error rather
than silently dropping detections. Input/provenance manifests contain no primer
observations. A blank primer can
mean no hit or no tracking row; `primer_tracking_record_found` distinguishes these.

Combined tables:

- `all_reads.tsv`: concatenated read tables, retaining `source_folder`, `source_name`,
  full `source_file`, `dataset`, and `PCR_plate` for every row.
- `all_species_summary.tsv`: concatenated species summaries with folder/file provenance.
- `barcode_analysis_summary.tsv`: per-folder read/assignment/primer counts and forward/reverse detected-name lists.
- `species_by_folder.tsv`, `contamination_summary_by_folder.tsv`: reporting tables.
- `input_manifest.tsv`, `taxon_plate_prevalence.tsv`: inputs and flag provenance.

The parser handles the large files one barcode at a time, with temporary intermediates
under `/tmp`; it does not keep all read tables in memory simultaneously.

## Export cleaned reads by taxon

```bash
python code/export_sintax_fasta.py \
  --taxon Insecta \
  --folders barcode01 barcode02 barcode15 \
  --analysis-dir data/sintax_analysis/mergedseq5_20261021 \
  --output data/insecta_barcode01_02_15.fasta

python code/export_sintax_fasta.py \
  --taxon Tipulidae \
  --folders barcode15 \
  --analysis-dir data/sintax_analysis/mergedseq5_20261021 \
  --output data/barcode15_tipulidae.fasta
```

Names are case-insensitive exact matches against any reported SINTAX taxon label
at any rank with confidence >=0.80 on the matched label.
A complete MIDORI label or numeric taxid is also accepted. Repeat `--taxon` to take
the union of several taxa. A name must actually appear in the SINTAX lineage;
for example, a superfamily omitted from the database cannot be inferred by this
script. For Tipuloidea, explicitly select its represented families with repeated
`--taxon Tipulidae --taxon Limoniidae --taxon Pediciidae --taxon Cylindrotomidae`.
Misspellings such as `Tipulidea` produce a no-match error and possible suggestions,
not a silent substitution.

**Selection and annotation use separate rules.** A read is selected if the queried
taxon name, complete label or taxid appears in its reported SINTAX lineage
with confidence >=0.80 on that label. Lower ranks may have lower confidence.
Labels with missing confidence or confidence below 0.80 do not qualify.
The FASTA header uses the **lowest trusted taxon**: the lowest rank with confidence
>=0.80. If no rank reaches 0.80, the header taxon is `unassigned`. The contamination
flag refers to this trusted header taxon, not necessarily the taxon used to query.

For example, querying `Tipulidae` selects a read with family confidence 0.90
and genus `Tipula` confidence 0.20. Its header uses `Tipulidae_<taxid>`.
Querying `Tipula` would not select that read.

 Selected query IDs are fetched exactly from the matching folder's **cleaned**
`filtered/asinput.fastq.gz`. The cleaned sequence orientation is preserved. Missing
selected IDs abort the export without publishing a partial FASTA. Existing FASTA
outputs are never overwritten. The input SINTAX file must match the analysed file's
path, size and modification time so flags cannot silently become stale.

Headers have six logical components:

```
folder-forwardprimer-reverseprimer-lowest_trusted_taxon-contamination_or_nocontamination-readID
```

Example:

```
>barcode15-bc1087-bc1100rc-Tipula_maxima_12345-nocontamination-<original-cleaned-read-ID>
```

Both tags are canonical primer names for that read's raw parent. `none` means a
tracking record exists but that direction has no detected tag; `unknown` means
there is no tracking record. The taxon is the **lowest rank with confidence >=0.80**,
not an unsupported species candidate. If no rank qualifies it is `unassigned`.
The exports require an analysis at cutoff 0.80. `taxapass/taxafail` remains available
in the read tables but is no longer included in FASTA headers.

Whitespace in labels becomes underscores; embedded dashes in folder/taxon components
are percent-escaped. Dashes inside the final read ID are preserved. Split a header
at most five times to recover the components. `build_phylo_read_fastas.py` forwards
the same exporter options for compatibility. No consensus or amino-acid filtering occurs.

## Classification, only when explicitly wanted

The classification entry point remains available independently:

```bash
bash code/run_metabarcoding_pipeline.sh --dry-run data/mergedseq5/barcode01
bash code/run_metabarcoding_pipeline.sh data/mergedseq5/barcode01
```

It reuses existing filtered combined reads, decompressed MIDORI FASTA, and UDB.
If the UDB is absent it builds one; if the FASTA is absent it decompresses the `.gz`
archive first. Existing classifier outputs are skipped. `config/pipeline.env`
contains the paths and settings. The existing classifier command uses cutoff 1;
post-classification analysis independently uses recorded confidence >=0.8. Setting
`SINTAX_CUTOFF=0.8` changes future classifier runs only. Nothing in the notebook
runner invokes classification, even when barcode outputs are missing.

## Validation

```bash
python code/test_read_analysis.py -v
python code/test_sintax_migration.py -v
```

Tests use temporary fixtures, including threshold boundaries, distinct PCR-plate
counts, confidence medians, primer joins, missing outputs, clean-read export and
existing-file protection.

## Export by forward-primer tag

Use `code/export_primer_fasta.py` to select raw-parent best forward-tag hits, trace
them into cleaned reads, and export FASTA with the same SINTAX taxonomy and
contamination headers. See [primer export documentation](export_primer_fasta.md)
for the bc1087/barcode15 example and optional `--include-unclassified` search.

Module 2 reference selection requires either `--download` for NCBI COI records
or `--bold` for taxon-name matching against compressed FASTA headers in
`data/bold/`. Both feed the same VSEARCH, optional translation-filtering, MAFFT or kc-align alignment,
phylogeny and Baltic plotting steps. See the
[Module 2 source options](post_processing.md#choose-the-reference-dataset).
