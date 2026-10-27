# Deadwood insect COI metabarcoding

The pipeline has two modules:

1. **Cleaning and taxonomic assignment**: `code/run_metabarcoding_pipeline.sh`
   cleans reads and runs SINTAX; `code/run_sintax_notebooks.py` generates read-level
   tables, primer summaries and contamination flags. See the
   [Module 1 documentation](docs/pipeline_from_fastplong_onwards.md).
2. **Phylogenetics**: `code/postprocess_sintax.py` exports a selected taxon,
   selects NCBI COI references (`--download`) or taxon-matched local BOLD records (`--bold`), filters references to 500–2,000 nt and builds cluster consensus sequences at 98% identity with VSEARCH (local barcode/primer/taxon groups; references globally), optionally filters six-frame translations,
   aligns with MAFFT FFT-NS-2 or MACSE (`--macse`), estimates an IQ-TREE phylogeny, and renders a midpoint-rooted
   Baltic divergence-tree figure with supports, trusted-taxon/ID tip labels and scale bar. Mitogenomes are excluded
   by default; `--include-mitogenomes` allows them. See the
   [Module 2 documentation](docs/post_processing.md).

Oxford Nanopore cleaned-read taxonomic analysis using VSEARCH SINTAX and MIDORI COI
references. The post-classification workflow parses all taxonomic ranks/confidences,
assigns the lowest rank with confidence >=0.8, joins primer-tracking evidence,
and flags taxa occurring in >=80% of analysed PCR-plate codes.

Per-barcode and combined tab-delimited tables preserve query IDs, parent-read IDs,
source folder/file, and PCR-plate provenance. The taxon FASTA exporter selects reads matching any reported SINTAX taxon label,
with confidence >=0.80 on the queried label. Exported cleaned sequences are annotated with the lowest
trusted taxon (lowest rank with confidence >=0.80), or `unassigned` if none qualifies.

```bash
conda activate taxonomy
python code/run_sintax_notebooks.py --output-dir data/sintax_analysis/my_new_run
```

This command executes post-classification notebooks only. Missing classifier outputs
are reported and excluded. Choose a new output directory for each run.

See [workflow documentation](docs/pipeline_from_fastplong_onwards.md) for definitions,
outputs, confidence/contamination rules, FASTA export and validation.

For taxon-specific COI reference retrieval, grouped/global VSEARCH 98% nucleotide clustering, six-frame
mitochondrial translation filtering, alignment and phylogeny, see the
[post-processing module](docs/post_processing.md).

Module 2 supports `--resume` in the same output directory: retained sequences
restart at the chosen aligner, completed alignments restart at IQ-TREE, and earlier saved
download/extraction/cluster outputs are reused when needed. Batched GenBank XML
is parsed locally without repeating completed metadata downloads.

To export cleaned reads selected by the best forward-primer tag on their raw
parent reads (for example bc1087 in barcode15), see the
[standalone primer exporter](docs/export_primer_fasta.md). Optional searching of
`unclassified` is off by default.

Local FASTA headers contain folder, forward primer, reverse primer, the lowest
taxon trusted at confidence >=0.80, contamination flag and cleaned read ID.
Use a new output directory when changing biological settings or grouping policy.

Taxon queries use reported labels with confidence >=0.80 on the queried rank; FASTA headers use the lowest
trusted taxon at confidence >=0.80. These are separate rules. Use a new Module 2
output directory for analyses using changed selection settings.

Module 2 requires exactly one reference-source flag: `--download` (NCBI Entrez)
or `--bold` (offline header matching in `data/bold/*.fasta.gz`). BOLD mode accepts
a taxon name and labels reference tips with the available taxon and BOLD process
ID. Both modes use the same clustering, translation filter and phylogeny steps;
see the [source options and examples](docs/post_processing.md#choose-the-reference-dataset).

For nuclear/noncoding loci, use default MAFFT with `--skip-translation` and
an appropriate marker-specific database. The NCBI download option still targets
COI. See [alignment and translation options](docs/post_processing.md).
