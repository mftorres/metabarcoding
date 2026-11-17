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
  --analysis-dir data/read_processing/sequencing2_run01/analysis \
  --output-dir data/postprocessing/Tipulidae_barcode15_bold_updated_run01 \
  --threads 6
```

Replace the taxon, barcode, completed analysis path and new results folder with
your own choices. If an old command explicitly names an older database after
`--bold-db`, change that argument to the new file; renaming another file does
not update the command automatically. The HTML command builder's **BOLD database
path** field should contain `data/bold/BOLD_filtered.fasta.gz`.

Without `--bold-db`, the script scans **all** supported files in `data/bold/`.
Keep only the intended prepared snapshot there when using that default; raw,
filtered and old copies together can include the same records more than once.
Specifying the prepared file explicitly is the clearest way to select a release.
Start a new Module 2 run when changing the database; do not resume a run using
a different snapshot.

Empty-record removal is a preparation step, not a species-identification check.
Module 2 still keeps only named species identifications from BOLD, filters
references to 500–2,000 nucleotides, and clusters them at 98% within each taxonomic
label. Local reads use their separate dominant-label selection and 97% clustering.


## BOLD headers and downstream selection

Preserve the full original FASTA headers. The expected structure is:

```text
>processID|marker|country|comma-separated lineage
```

Use a taxon name, not a numeric taxid. Header matching includes all matching
markers, so use a marker-specific snapshot if needed. The pipeline generates
selected references and metadata inside each run directory; those outputs do not
belong in this database folder. See the [Postprocessing guide](../../docs/post_processing.md).

Postprocessing extracts BOLD records as described above, then clusters at 98%
identity within each exact taxonomic suffix after `__taxon__` in the prepared
FASTA headers. Different species labels never share a consensus, even if their
sequences are identical, so their positions remain visible in the tree. See the
[Postprocessing guide](../../docs/post_processing.md) for the separate dominant-label
selection of local reads and the extraction audit columns.

BOLD references must have a named species identification. In full BOLD lineages,
the eighth comma-separated field is species (empty positions are preserved);
the optional ninth subspecies field is not used as the species label. Shortened
lineages are accepted only when their final nonempty entry is a named binomial.
Genus-only, missing and provisional identifications (`sp.`, `spp.`, `cf.`, `aff.`)
are excluded. Named binomials, optional subgenus annotations and trinomial names
are accepted. If no qualifying records match, extraction stops with an error.

## Review the chosen reference set

BOLD is the preferred starting source for this COI workflow. Broad NCBI downloads
can produce large, heterogeneous sets requiring more record-level review; this
does not mean all NCBI records are uncurated or all BOLD identifications are
correct. Which BOLD species and records to keep is the user’s scientific decision.
Review record quality and the trimmed alignment, document exclusions, and
re-estimate the tree if the alignment or sequence set changes. See the
[alignment and tree review instructions](../../docs/post_processing.md#review-the-trimmed-alignment-before-interpreting-the-tree).
