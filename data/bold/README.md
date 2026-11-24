# BOLD reference database for Postprocessing

Run commands from the repository root with the `taxonomy` environment active.
The database stays local and is excluded from Git; commit this README instead.
Postprocessing needs it only for `--bold`; `--download` obtains references from NCBI.

## Download and prepare an updated BOLD database

Do this before running Module 2 with `--bold`. A database downloaded earlier does
not update itself. The following steps prepare a new snapshot without editing
the pipeline code.

### 1. Download the sequence file

Open the [latest BOLD data package](https://bench.boldsystems.org/index.php/datapackages/Latest)
in your web browser. If that address is unavailable, open the
[BOLD Data Packages page](https://boldsystems.org/data/data-packages/) and choose
**Latest**. Download the **FASTA sequence file**, rather than the TSV data table
or JSON metadata. If BOLD requests sign-in, follow its access instructions.

Save the download in your normal Downloads folder and note its release date,
original filename and source link. “Latest” changes over time; keeping these
notes lets you identify which release you used. If the download is an outer ZIP
or TAR archive, extract it and locate the sequence FASTA first. An actual
`.fasta.gz` sequence file can stay compressed. Renaming a TSV file to `.fasta.gz`
does not convert it into sequence data.

### 2. Open Ubuntu and locate the download

Go to your metabarcoding project folder and activate the environment. Change
the first line if the project is stored elsewhere:

```bash
cd ~/metabarcoding
conda activate taxonomy
seqkit version
```

If you downloaded through Windows, list the Downloads folder from Ubuntu:

```bash
ls -lh "/mnt/c/Users/USERNAME/Downloads/"
```

Replace `USERNAME` with the Windows account folder name. For a download made
inside Linux, use `ls -lh ~/Downloads/` instead. Find the actual FASTA filename.
In the next command, **replace the entire example input path** with that file.
`DOWNLOAD_DATE` below is a placeholder, not a required part of the name.

### 3. Remove records that have a header but no DNA sequence

A FASTA record starts with a line beginning with `>`, followed by its DNA letters.
Some BOLD records have a header but no letters. Use SeqKit to remove those empty
records before the pipeline reads the database:

```bash
mkdir -p data/bold data/bold_staging
seqkit seq --min-len 1 \
  "/mnt/c/Users/USERNAME/Downloads/BOLD_Public.DOWNLOAD_DATE.fasta.gz" \
  -o data/bold_staging/BOLD_download.filtered.fasta.gz
```

The input can be plain FASTA or gzipped FASTA. SeqKit writes compressed output
because the output name ends in `.gz`; you do not need to decompress the input
or compress the result separately. Keep input and output as different files.
The staging folder is outside `data/bold/`, so the pipeline will not accidentally
scan a file while you are preparing it.

`--min-len 1` means “keep sequences at least one DNA letter long.” It removes
zero-length records, keeps the full headers of retained records, and does not
select species, remove duplicates or apply the later reference-length filter.
See the [SeqKit command documentation](https://bioinf.shenwei.me/seqkit/usage/#seq).

Check the filtered file before continuing:

```bash
seqkit stats data/bold_staging/BOLD_download.filtered.fasta.gz
gzip -t data/bold_staging/BOLD_download.filtered.fasta.gz
```

The statistics should show at least one sequence and a minimum sequence length
of at least 1. `gzip -t` normally prints nothing when successful; an error means
you should fix or repeat preparation before proceeding. For large databases,
filtering and these checks can take time. Keep the terminal open and wait for
the prompt to return. If SeqKit reports an error, do not continue with its output.

### 4. Give the prepared database a clear name and put it in the right folder

Use this command to move and rename the filtered file:

```bash
mv -i data/bold_staging/BOLD_download.filtered.fasta.gz data/bold/BOLD_filtered.fasta.gz
ls -lh data/bold/BOLD_filtered.fasta.gz
```

The file is now called `BOLD_filtered.fasta.gz` inside `data/bold/`.
If that filename already exists, `mv -i` asks before replacing it. Enter `n` to
keep the existing file, or `y` only if you intend to replace it with the new one.
If you chose `n`, the new file remains in the staging folder. Keep a dated copy
or move the old version outside `data/bold/` if you need to preserve it.

**There is no required date or fixed database basename in the scripts.**
`BOLD_filtered.fasta.gz` is the standard name used in these instructions, not a
hard-coded requirement. The default folder search accepts `.fasta.gz`, `.fa.gz`
and `.fna.gz`. You can keep the real release date in your filename instead and
use that name in `--bold-db`. Do not rename a new release with an old release's
date simply to match an older example.

### 5. Tell Module 2 to use the prepared file

Add `--bold-db data/bold/BOLD_filtered.fasta.gz` to your BOLD command. For example:

```bash
python code/postprocess_sintax.py --bold \
  --bold-db data/bold/BOLD_filtered.fasta.gz \
  --taxon Tipulidae --folders barcode15 \
  --analysis-dir data/read_processing/DATASET_run01/analysis \
  --output-dir data/postprocessing/Tipulidae_barcode15_bold_updated_run01 \
  --threads 6
```

Replace the taxon, barcode, completed analysis path and new results folder with
your own choices. If an old command explicitly names an older database after
`--bold-db`, change that argument to the new file; renaming another file does
not update the command automatically. In the HTML guide's command builder, enter the same path in the **BOLD database** field.

Without `--bold-db`, the script scans **all** supported files in `data/bold/`.
Keep only the intended prepared snapshot there when using that default; raw,
filtered and old copies together can include the same records more than once.
Specifying the prepared file explicitly is the clearest way to select a release.
Start a new Module 2 run when changing the database; do not resume a run using
a different snapshot.

Empty-record removal is a preparation step, not a species-identification check;
see [How Module 2 selects BOLD records](#how-module-2-selects-bold-records).


## How Module 2 selects BOLD records

Keep the full original FASTA headers. Their structure is:

```text
>processID|marker|country|comma-separated lineage
```

Module 2 (`--bold --taxon NAME`) selects every record whose header contains `NAME`
(not case-sensitive), anywhere in the header. `--taxon Tipula` therefore also selects
*Tipulidae* records; this is intended, so that a genus search brings in its family as
context. Use a taxon name, not a numeric taxid.

Only records with a **named species** are kept: the eighth lineage field (species)
must be a proper binomial (optionally with subgenus or subspecies). Genus-only,
empty and provisional names (`sp.`, `spp.`, `cf.`, `aff.`, `nr.` …) are excluded.
Selected records are then length-filtered (500–2,000 nt) and clustered at 98%
identity within each species; different species never share a consensus. The
selected references and their metadata are written into each run folder, not here.

## Other markers in the BOLD file

The BOLD file is not COI-only. In the 16-Oct-2026 snapshot, 20.5 million records are
COI-5P (the barcode region) and about 1.5 million (6.5%) are other markers, for example
ITS, COI-3P, 18S, matK, rbcL, CYTB, 16S, 28S, COII and the ND genes. Module 2 does not
filter by marker, so records of other markers whose header contains your taxon name
are selected as references.

- The protein check removes most **non-coding** markers (ITS, 16S, 18S, 28S), because
  they rarely contain 200 amino acids without a stop codon.
- **Protein-coding** markers (COI-3P, CYTB, COII, ND genes) can pass it. kc-align with a
  COI guide normally rejects them, because they do not align to COI.
- **If you use MAFFT, or `--skip-translation`, other markers can remain in the
  alignment and the tree.** Look for references that do not line up with the COI
  sequences during review, and check their marker in `bold_reference_metadata.tsv`
  (column `marker`).

To exclude other markers from the start, you can prepare a COI-5P-only file:

```bash
seqkit grep --by-name --use-regexp --pattern '\|COI-5P\|' \
  data/bold/BOLD_filtered.fasta.gz -o data/bold_staging/BOLD_COI-5P.fasta.gz
```

and use it with `--bold-db data/bold_staging/BOLD_COI-5P.fasta.gz` (or move it into
`data/bold/` instead of the full file).

## Review the chosen reference set

BOLD is the preferred starting source for this COI workflow. Broad NCBI downloads
can produce large, heterogeneous sets requiring more record-level review; this
does not mean all NCBI records are uncurated or all BOLD identifications are
correct. Which BOLD species and records to keep is the user’s scientific decision.
Review record quality and the trimmed alignment, document exclusions, and
re-estimate the tree if the alignment or sequence set changes. See the
[review instructions](../../docs/post_processing.md#c-review-the-results).
