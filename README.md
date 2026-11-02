# Deadwood insect COI metabarcoding

Open the [offline browser guide](docs/metabarcoding-guide.html) for WSL/Conda setup,
copyable commands, a left-hand process menu and a Postprocessing command builder.
Double-click the HTML file; no hosting, compilation or server is required.
The [taxonomy environment recipe](config/taxonomy-environment.yml) records current
direct dependency versions for Linux x86_64 / WSL2.

Get the project from [GitHub](https://github.com/mftorres/metabarcoding):

```bash
cd /mnt/c/Users/USERNAME/Documents
git clone --branch devfull https://github.com/mftorres/metabarcoding.git
cd metabarcoding
```

Replace `USERNAME` with your Windows account folder name. See the
[Ubuntu/WSL setup instructions](docs/setup.md) for the Linux Miniconda download,
installation, Windows Downloads/Documents access, and environment creation.

For faster file I/O under WSL, consider running analyses in the local Linux
filesystem (for example `~/metabarcoding`) rather than `/mnt/c/Users/USERNAME/`.
Keep inputs, databases and outputs there; choose the location before starting
analyses. See the [setup performance note](docs/setup.md#performance-run-analyses-on-the-linux-filesystem).

Before running, download the external reference databases. The folder READMEs
provide source links, release selection, exact destination paths, compression
instructions and BOLD empty-record cleanup with SeqKit:

- [MIDORI SINTAX setup (Read processing)](data/sintax/README.md)
- [BOLD FASTA setup (Postprocessing `--bold`)](data/bold/README.md)

Database files remain local; the setup READMEs are version controlled.

Run a dataset by passing its parent folder once (a bare `sequencing2` resolves to
`data/sequencing2`; any dataset name is accepted). Barcode folders and `unclassified` are discovered automatically:

```bash
# If reads are nested in pass/fail folders, merge them first; sources are preserved.
python code/merge_pass_fail.py sequencing2 --dry-run
python code/merge_pass_fail.py sequencing2
python code/read_processing.py sequencing2 --output-dir data/read_processing/sequencing2_run01
python code/postprocess_sintax.py sequencing2 --bold --taxon Tipulidae \
  --analysis-dir data/read_processing/sequencing2_run01/analysis \
  --output-dir data/postprocessing/Tipulidae_run01 --threads 8
```

Classification outputs are separated by dataset under `data/sintax/DATASET/`.
`unclassified` is included in read analyses but excluded from PCR-plate prevalence.
See the [full workflow](docs/pipeline_from_fastplong_onwards.md#run-both-modules-using-the-dataset-parent)
for prerequisites, path overrides and output handling.

The pipeline has two modules:

1. **Read processing**: `code/read_processing.py DATASET` runs
   cleaning, SINTAX, primer tracking and all reports. The lower-level
   `code/run_metabarcoding_pipeline.sh` and `code/run_sintax_notebooks.py` remain
   available for individual stages. See the
   [Read processing documentation](docs/pipeline_from_fastplong_onwards.md).
2. **Postprocessing**: `code/postprocess_sintax.py` exports a selected taxon,
   selects NCBI COI references (`--download`) or taxon-matched local BOLD records (`--bold`), filters references to 500–2,000 nt and builds cluster consensus sequences at 98% identity with VSEARCH (local barcode/primer/taxon groups; references globally), optionally filters six-frame translations,
   aligns with MAFFT FFT-NS-2 or kc-align (`--kcalign --kcalign-reference PATH`), trims alignment columns with ClipKIT smart-gap, estimates an IQ-TREE maximum-likelihood phylogeny without bootstrap replicates, and renders a midpoint-rooted
   Baltic divergence-tree figure with supports, trusted-taxon/ID tip labels and scale bar. Mitogenomes are excluded
   by default; `--include-mitogenomes` allows them. See the
   [Postprocessing documentation](docs/post_processing.md).

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

Postprocessing supports `--resume` in the same output directory: retained sequences
restart at the chosen aligner, completed alignments restart at ClipKIT (or IQ-TREE when trimming is complete), and earlier saved
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
trusted taxon at confidence >=0.80. These are separate rules. Use a new Postprocessing
output directory for analyses using changed selection settings.

Postprocessing requires exactly one reference-source flag: `--download` (NCBI Entrez)
or `--bold` (offline header matching in `data/bold/*.fasta.gz`). BOLD mode accepts
a taxon name and labels reference tips with the available taxon and BOLD process
ID. Both modes use the same clustering, translation filter and phylogeny steps;
see the [source options and examples](docs/post_processing.md#choose-the-reference-dataset).

For nuclear/noncoding loci, use default MAFFT with `--skip-translation` and
an appropriate marker-specific database. The NCBI download option still targets
COI. See [alignment and translation options](docs/post_processing.md).

The `code/` folder contains the two modules, their runtime helpers and source
notebooks. Standalone sequence extraction remains available through
`code/export_sintax_fasta.py` (taxon) and `code/export_primer_fasta.py`
(barcode/forward primer); `code/export_sintax_taxa.py` is the taxon entry point
used internally by Postprocessing. Saved executed notebooks are analysis outputs.
