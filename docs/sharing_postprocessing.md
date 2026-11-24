# Minimum inputs to share for Postprocessing

The tar-ready lists beside this guide were generated from the completed analysis
`data/sintax_analysis/mergedseq5_20261021/`. They include every analysed folder
(barcode01–barcode21), without restricting the local data to any particular taxon.
They support any taxon represented in those SINTAX assignments at the required
confidence, subject to Postprocessing's sequence/reference requirements.

This analysis does not include barcode22–barcode24 or `unclassified`. To share
those folders for Postprocessing, first include them in a completed Read
processing analysis and build an equivalent input list for that analysis.
Adding raw files alone does not provide the required assignments and annotations.

## Files to upload

Use `docs/postprocessing-share-data.list` as the primary file list. Its 68 paths
are relative to the project root, one file per line:

| Required files | Why Postprocessing needs them | Current size before tar compression |
| --- | --- | ---: |
| 21 × `data/mergedseq5/BARCODE/filtered/asinput.fastq.gz` | Actual cleaned sequences from which any requested taxon can be extracted. | 1.924 GiB |
| 21 × `data/sintax/BARCODE/BARCODE_sintaxSP.tsv` | Per-query taxonomy and confidences used to select taxa. | 2.036 GiB |
| 21 × `data/sintax_analysis/mergedseq5_20261021/BARCODE/BARCODE_reads.tsv` | Saved primer annotations for extraction and taxon/read labels for tree plotting. | 4.433 GiB |
| `analysis.complete` in that analysis folder | Proof of completed table generation required by the pipeline. | Small |
| `analysis_parameters.json` in that analysis folder | Dataset/classification paths, analysed folders and confidence settings. | Small |
| `input_manifest.tsv` in that analysis folder | Original classifier paths, sizes and timestamps used to validate inputs. | Small |
| `primer_name_mapping.json` in that analysis folder | Canonical primer-name mapping. | Small |
| `taxon_plate_prevalence.tsv` in that analysis folder | Existing prevalence contamination annotations. | Small |

The total is **9,013,334,678 bytes (about 8.394 GiB)** before tar compression.
The TSVs will be compressed by tar/gzip; the FASTQs are already compressed, so
expect little additional compression for that part. Actual archive size is not
known until compression finishes.

The per-barcode read tables contain extra reporting columns, but these are the
existing files consumed by the current pipeline. They also avoid needing to ship
the much larger original primer-hit tables. Do not omit them: plotting uses them
in addition to extraction's primer lookup. A specially prepared smaller copy can
retain only `readid`, `raw_ont_read_id`, `forward_PCR_primer`, `reverse_PCR_primer`,
`primer_tracking_record_found`, `accepted_species` and `assigned_name`, while
keeping every row. The supplied list uses the original tables without rewriting
analysis data.

## Code, environment and reference-source choices

Collaborators can clone the code rather than upload another copy:

```bash
git clone --branch devfull https://github.com/mftorres/metabarcoding.git
cd metabarcoding
```

Install the `taxonomy` environment using `config/taxonomy-environment.yml` and
follow [setup](setup.md). If GitHub access is unavailable,
`docs/postprocessing-share-code.list` lists the runtime code and environment
recipe to add to the same archive. It includes the kc-align helpers as well as
the default MAFFT path.

Reference inputs depend on the chosen mode:

| Mode | Additional requirement | Upload list |
| --- | --- | --- |
| `--download` with default MAFFT | Internet access plus Entrez Direct; references are fetched during the run. No local comparison database is required. | No additional data list. |
| `--bold` | The full BOLD snapshot, or a compatible marker-specific snapshot containing every taxon the collaborator needs. | `docs/postprocessing-share-bold.list` contains your local `data/bold/BOLD_Public.16-Oct-2026.fasta.gz` (about 2.25 GiB). Downloading it separately avoids including it in each upload. |
| `--kcalign` | A suitable single-sequence coding reference supplied via `--kcalign-reference`. | `docs/postprocessing-share-kcalign.list` lists the three small guides currently in `data/postprocessing/coi_refs/`. These are optional and are not universal references for every taxon/marker. |

See the [BOLD setup guide](../data/bold/README.md) for downloads. A taxon-specific
reference subset only supports the taxa represented in that subset; it is not a
replacement for a broad database when collaborators need arbitrary taxa.

## Create the tar.gz

Run from the project root. These commands package files; they do not rerun any
analysis. Use PAX format so the classifier nanosecond modification timestamps
needed by the manifest checks are preserved.

```bash
tar --format=pax -czf postprocessing-inputs.tar.gz \
  --verbatim-files-from --files-from=docs/postprocessing-share-data.list
```

To include runtime code and the available kc-align guides as well:

```bash
tar --format=pax -czf postprocessing-inputs-with-code.tar.gz \
  --verbatim-files-from \
  --files-from=docs/postprocessing-share-data.list \
  --files-from=docs/postprocessing-share-code.list \
  --files-from=docs/postprocessing-share-kcalign.list \
  docs/sharing_postprocessing.md
```

Add `--files-from=docs/postprocessing-share-bold.list` only if you intend to ship
the BOLD database in this archive. Do not recompress the BOLD download first;
it is already gzipped.

For a per-file upload limit, split the compressed archive. For example, these
parts stay below a 1 GB limit (choose a smaller `-b` value for a lower limit):

```bash
split -b 900M -d -a 3 postprocessing-inputs.tar.gz postprocessing-inputs.tar.gz.part-
```

The recipient needs **all** parts and can reconstruct the archive with:

```bash
cat postprocessing-inputs.tar.gz.part-* > postprocessing-inputs.tar.gz
```

For a complete illustrated-by-commands recipient route, including folder creation,
extraction and the one-time path update, see the [Postprocessing archive walkthrough](post_processing.md#if-you-downloaded-postprocessing-inputstargz-start-here)
or the [self-contained HTML instructions](metabarcoding-guide.html#shared-inputs).

## Relocate the saved input paths on the recipient's computer

Extract from inside the cloned/project folder, preferably on the local Linux
filesystem. The paths stored in your existing analysis point to the sender's
computer. Even if your project folder is also named `metabarcoding`, its absolute
path will usually differ. Extracting the archive does not update those saved
paths; run the verified relocation command below before Postprocessing.

```bash
tar -xzf /path/to/postprocessing-inputs.tar.gz
```

The following one-time command is for the specific dataset/analysis in these
lists. It verifies classifier sizes and original timestamps, checks required read
annotation columns, backs up the two original metadata files, and then updates
only input paths. It preserves confidence, taxonomy, primer mappings, prevalence
and all recorded size/timestamp checks. It does not change the sequence data.

```bash
python - <<'PY'
import csv
import json
import shutil
from pathlib import Path

root = Path.cwd().resolve()
analysis = root / 'data/sintax_analysis/mergedseq5_20261021'
dataset = root / 'data/mergedseq5'
sintax = root / 'data/sintax'
params_path = analysis / 'analysis_parameters.json'
manifest_path = analysis / 'input_manifest.tsv'
params = json.loads(params_path.read_text())
with manifest_path.open() as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    columns = reader.fieldnames
    rows = list(reader)
required = {'readid', 'raw_ont_read_id', 'forward_PCR_primer',
            'reverse_PCR_primer', 'primer_tracking_record_found',
            'accepted_species', 'assigned_name'}
selected = set(params['analysed_folders'])
assert selected.issubset({row['source_folder'] for row in rows})
for row in rows:
    folder = row['source_folder']
    if folder not in selected:
        continue
    source = sintax / folder / f'{folder}_sintaxSP.tsv'
    clean = dataset / folder / 'filtered/asinput.fastq.gz'
    table = analysis / folder / f'{folder}_reads.tsv'
    assert clean.is_file(), clean
    stat = source.stat()
    assert stat.st_size == int(row['sintax_bytes']), source
    assert stat.st_mtime_ns == int(row['sintax_mtime_ns']), (
        f'{source}: timestamp changed; extract the original PAX archive '
        'on Linux rather than replacing the saved validation values')
    with table.open() as handle:
        fields = next(csv.reader(handle, delimiter='\t'))
    assert required.issubset(fields), table
    row['sintax_file'] = str(source.resolve())
    row['clean_reads_file'] = str(clean.resolve())
for path in (params_path, manifest_path):
    backup = path.with_name(path.name + '.before_relocation')
    if not backup.exists():
        shutil.copy2(path, backup)
params['merged_dir'] = str(dataset)
params['sintax_root'] = str(sintax)
params_path.write_text(json.dumps(params, indent=2) + '\n')
with manifest_path.open('w', newline='') as handle:
    writer = csv.DictWriter(handle, columns, delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)
print('Validated and relocated Postprocessing inputs.')
PY
```

The old `primer_hits` path remains as provenance in the parameters. Older
analyses may also record a retired PCR-map path; it is not read by Postprocessing.
Postprocessing does not need those external files when the complete per-barcode
read tables and saved primer mapping/prevalence are provided as listed here.
This compact package is for Postprocessing, not for rerunning raw-read primer
tracking or the standalone raw-parent primer exporter.

Now a collaborator can select any included folder/taxon, for example:

```bash
conda activate taxonomy
python code/postprocess_sintax.py mergedseq5 --download --taxon Tipulidae \
  --analysis-dir data/sintax_analysis/mergedseq5_20261021 \
  --output-dir data/postprocessing/Tipulidae_collaborator_run01 --threads 8
```

Use `--folders barcode01 barcode02` instead of the dataset argument for a subset.
Use `--bold` and the prepared BOLD snapshot if that is the desired reference source.
Use a new output directory for a new run.

## Files that are not needed in this upload

Do not include raw barcode FASTQs, per-file filtered FASTQs when `asinput.fastq.gz`
is present, the original classification database, raw SeqKit hits or the original primer best-hit table,
`all_reads.tsv` (duplicates the per-barcode tables), species/primer summary tables,
figures or executed notebooks. Existing `data/postprocessing/*run*/` outputs,
particularly IQ-TREE distance matrices, model caches and temporary archives, are
not inputs for a new arbitrary-taxon run. The optional coding guides under
`coi_refs/` are the exception for kc-align.

## Barcode-only analyses and legacy sharing lists

New Read processing analyses identify samples by barcode folder and write
`taxon_barcode_prevalence.tsv`; they need no PCR/pool mapping file. The explicit
file list in this guide still names the existing 2026-10-21 analysis's
`taxon_plate_prevalence.tsv`. That is intentional: Postprocessing supports that
unchanged legacy file and its saved flags without requiring the retired mapping.
When preparing a list for a new analysis, include its
`taxon_barcode_prevalence.tsv` instead. Do not rename or regenerate older results
just to use the current code.


## Copy-and-send instructions for collaborators

Follow **Route A only** below. It supplies the processed inputs; do not run
Module 1 or Route B afterwards. These instructions focus on BOLD and kc-align for COI.


**Use only this route if you received the input archive.** It already supplies
processed inputs for the included barcode01–barcode21 samples. You still need the
code, software environment, BOLD database and a prepared COI coding guide.
**Do not run Module 1. Do not follow Route B afterwards.**

### A.1 Choose the Linux project location

All commands run in **Ubuntu**, not PowerShell. A path is a folder address.
`~` means your Linux home folder. Prefer **`~/metabarcoding`** for the project and
its data: Linux storage usually handles the many analysis files faster than the
Windows drive mounted under `/mnt/c/`. You can keep downloaded archives in Windows
Downloads and extract them directly into Linux.

To enter Windows Downloads from Ubuntu, replace `USERNAME` with your Windows
account folder name:

```bash
cd "/mnt/c/Users/USERNAME/Downloads"
ls -lh
```

To return to the Linux project, run:

```bash
cd ~/metabarcoding
pwd
```

`cd` changes where subsequent commands run; it does **not** copy or move files.
`pwd` prints your current folder. To view a Linux folder in Windows File Explorer,
run `explorer.exe .` from Ubuntu, or browse through `\\wsl$`, your Ubuntu
installation, `home`, your Linux username and `metabarcoding`.

For a new project, clone directly into your Linux home:

```bash
cd ~
git clone --branch devfull https://github.com/mftorres/metabarcoding.git
cd ~/metabarcoding
ls code config
```

Git creates `~/metabarcoding`. If this repository is already there, skip cloning
and use `cd ~/metabarcoding`; do not clone again inside it. If that path contains
something else, choose another short location and use that location consistently.

If you have a project in Windows and want a Linux copy, do this **only when
`~/metabarcoding` does not already exist**, before running your own Module 1:

```bash
cp -a "/mnt/c/Users/USERNAME/Documents/metabarcoding" ~/metabarcoding
cd ~/metabarcoding
```

Replace the Windows path with the real one. `cp -a` copies the folder and preserves
its file attributes; it leaves the Windows original in place. Choose either
cloning or copying, not both. Moving an already completed analysis also requires
its saved absolute paths to be updated; simply using `cd` or copying does not fix
those paths. The Route A relocation step below is specific to the shared archive.

### A.2 Activate the software

From the project root, create the environment **only if it has not been
installed**:

```bash
conda env create --name taxonomy --file config/taxonomy-environment.yml
conda activate taxonomy
```

If it already exists, run only `conda activate taxonomy`. Activate it whenever you
open a new terminal for analysis. If Git or Conda is not available, finish the
[computer setup instructions](setup.md) before continuing.

### A.3 Extract the archive into the project root

Save `postprocessing-inputs.tar.gz` in Windows Downloads or your Linux Downloads.
The archive is roughly 2.5 GB in this checkout and expands to about 9 GB; allow
additional space for BOLD, software and results. Check space with `df -h ~`.
If the sender split it into parts, download every part and join them using the
[sharing instructions](sharing_postprocessing.md) before extracting.

From **`~/metabarcoding`**, run:

```bash
cd ~/metabarcoding
mkdir -p data
tar --keep-old-files -xzf "/mnt/c/Users/USERNAME/Downloads/postprocessing-inputs.tar.gz"
```

Replace `USERNAME`. For Linux Downloads, use `~/Downloads/postprocessing-inputs.tar.gz`
instead. This decompresses and extracts in one step and leaves the downloaded
archive in place. The archive already contains paths starting with `data/`.
**Do not extract from inside `data/`**: that would create `data/data/...`.
Nested sample folders are created automatically.

`--keep-old-files` refuses to overwrite existing files. If extraction reports an
error or a collision, stop; some files may already have been extracted. Use a
fresh project location instead of analysing a mixture of old and new inputs.
These instructions are for the **data-only** archive, not one that also includes code.

Check the extracted structure:

```bash
ls data/mergedseq5/barcode01/filtered/asinput.fastq.gz
ls data/sintax/barcode01/barcode01_sintaxSP.tsv
ls data/sintax_analysis/mergedseq5_20261021/analysis.complete
```

The first file is cleaned DNA reads; the second is saved classification output;
the third marks the completed analysis. Keep all the extracted inputs. Module 2
does not open the original classification database, but it still needs these
saved classification results. If your archive has different dataset names, use
its accompanying paths rather than guessing replacements.

### A.4 Update the saved paths for this computer

The analysis records paths from the sender's computer. Extraction alone does not
update them. From `~/metabarcoding` with `taxonomy` active, copy this **entire block**,
including the last `PY` line. It checks sizes, timestamps and read annotations,
backs up metadata, and updates only the saved input paths:

```bash
python - <<'PY'
import csv
import json
import shutil
from pathlib import Path

root = Path.cwd().resolve()
assert (root / 'code/postprocess_sintax.py').is_file(), 'First cd into the project folder'
analysis = root / 'data/sintax_analysis/mergedseq5_20261021'
dataset = root / 'data/mergedseq5'
sintax = root / 'data/sintax'
params_path = analysis / 'analysis_parameters.json'
manifest_path = analysis / 'input_manifest.tsv'
params = json.loads(params_path.read_text())
with manifest_path.open() as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    columns = reader.fieldnames
    rows = list(reader)
required = {'readid', 'raw_ont_read_id', 'forward_PCR_primer',
            'reverse_PCR_primer', 'primer_tracking_record_found',
            'accepted_species', 'assigned_name'}
selected = set(params['analysed_folders'])
assert selected.issubset({row['source_folder'] for row in rows})
for row in rows:
    folder = row['source_folder']
    if folder not in selected:
        continue
    source = sintax / folder / f'{folder}_sintaxSP.tsv'
    clean = dataset / folder / 'filtered/asinput.fastq.gz'
    table = analysis / folder / f'{folder}_reads.tsv'
    assert clean.is_file(), clean
    stat = source.stat()
    assert stat.st_size == int(row['sintax_bytes']), source
    assert stat.st_mtime_ns == int(row['sintax_mtime_ns']), (
        f'{source}: timestamp changed; extract the original PAX archive '
        'on Linux rather than replacing the saved validation values')
    with table.open() as handle:
        fields = next(csv.reader(handle, delimiter='\t'))
    assert required.issubset(fields), table
    row['sintax_file'] = str(source.resolve())
    row['clean_reads_file'] = str(clean.resolve())
for path in (params_path, manifest_path):
    backup = path.with_name(path.name + '.before_relocation')
    if not backup.exists():
        shutil.copy2(path, backup)
params['merged_dir'] = str(dataset)
params['sintax_root'] = str(sintax)
params_path.write_text(json.dumps(params, indent=2) + '\n')
with manifest_path.open('w', newline='') as handle:
    writer = csv.DictWriter(handle, columns, delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)
print('Validated and relocated Postprocessing inputs.')
PY
```

Wait for `Validated and relocated Postprocessing inputs.` If it fails, stop and
send the error to the maintainer. Do not edit validation values or timestamps to
bypass the check. Extracting on Linux preserves the precision these checks need.
This relocation is for the dataset in the supplied archive, not for arbitrary
Module 1 analyses. No biological analysis is repeated here.

### A.5 Download and prepare BOLD

BOLD is the preferred starting source for this COI workflow. It focuses on
barcode records with specimen information; broad NCBI downloads can produce
large, heterogeneous sets needing more record-level review. Neither source
guarantees correct identifications. Curated NCBI records are still useful,
including for the separate coding guide. The BOLD database and the single
kc-align coding guide are **two different inputs**.



Do this before running Module 2 with `--bold`. A database downloaded earlier does
not update itself. The following steps prepare a new snapshot without editing
the pipeline code.

#### Download the sequence file

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

#### Open Ubuntu and locate the download

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

#### Remove records that have a header but no DNA sequence

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

#### Give the prepared database a clear name and put it in the right folder

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


Use `data/bold/BOLD_filtered.fasta.gz` in the command in A.7. The explicit
`--bold-db` prevents the pipeline from scanning other snapshots in the folder.
Empty-record filtering does not identify species; Module 2 later selects named
species records and applies its 500–2,000 nucleotide reference filter.

### A.6 Choose and prepare the COI coding guide for kc-align

This is not automatically supplied by the archive or by the BOLD database.
Prepare it before running A.7. Use a trustworthy sequence from NCBI or BOLD,
verify the precise COI coding interval, exclude suspected NUMTs and remove
primers. These full instructions are repeated here so you can stay in Route A.

#### Choose a trustworthy source


Use an identified specimen from the same species, if available, or a closely
related species in the same genus or family that covers the same COI region.
Prefer a curated or published record with specimen/voucher information and
supporting evidence that it is mitochondrial COI. Do not automatically use the
first search result, an unverified sequence, or a local consensus you are trying
to identify. A file already in `coi_refs/` is not proof of suitability for your taxon.

**From NCBI:** search [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore/) for
your group and COI, COX1 or cytochrome c oxidase subunit I. Inspect the GenBank
record, its mitochondrial source, publication, specimen and **CDS** annotation
(the protein-coding DNA feature). Check its strand, exact coordinates,
`codon_start`, translation and genetic code. Download the nucleotide sequence of
that COI feature or a verified target-region extract, not the protein FASTA or an
entire mitochondrial genome containing other genes. An annotated mitochondrial
genome can be a source for extracting COI even though whole genomes are normally
excluded from Module 2's comparison-reference download. See
[NCBI's explanation of GenBank features](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/).

**From BOLD:** select a named-species COI record with credible specimen
identification, voucher/institution information and, where available, trace data
and a supporting publication or linked GenBank record. Check the marker and
sequence quality, and download the nucleotide FASTA. The bulk database header
alone does not establish the exact coding boundaries or rule out a NUMT. If the
BOLD barcode lacks coding annotation, align it to a trustworthy annotated
mitochondrial COI sequence to establish its frame and precise target interval,
or ask the maintainer for a prepared guide.

#### Exclude suspected NUMTs

A **NUMT** is a copy of mitochondrial DNA present in the nuclear genome. The guide
should be genuine mitochondrial COI, not a NUMT, pseudogene, contaminant or chimera.
Check the source annotation and evidence of mitochondrial origin. Compare its
protein translation and nucleotide alignment with independently trusted COI
records. Reject suspected NUMTs, unexpected frameshifts, internal stops or
unexplained insertions/deletions.

**No internal stops is necessary but not sufficient.** Some NUMTs retain an
apparently intact reading frame. A database name or a successful pipeline check
cannot certify that the guide is NUMT-free. If the evidence is uncertain, choose a
better documented record rather than calling it verified. This limitation is
supported by [research on insect COI NUMTs](https://pmc.ncbi.nlm.nih.gov/articles/PMC10249859/).

#### Prepare the exact COI coding region, without primers

**Use the exact expected COI open reading frame (ORF) or the precisely defined
coding interval being analysed.** Establish its boundaries from annotation and
known codon positions, not by guessing a length. There is no single universal COI
length for all taxa and all assays.

- If your analysis uses the **complete COI ORF**, use the complete annotated COI
  coding sequence of the expected length for that record/taxon, in the correct
  orientation. A short COI-5P barcode is not a complete COI ORF.
- If your reads cover only a **COI barcode fragment**, use a verified, in-frame
  extract corresponding to that same region. It must start and finish on complete
  codon boundaries, with its exact expected coding length documented. Do not use
  an unrelated fragment or assume a fixed barcode length establishes the frame.
- Remove laboratory primers, primer tails, tags, adapters and flanking non-COI
  sequence. Use the known primer positions and coding annotation to define the
  interval; do not remove genuine COI bases solely because a short motif resembles
  a primer. The guide must contain only the intended biological coding region.
- Save **one nucleotide sequence** in a plain FASTA file: one `>` header followed
  by ungapped A/C/G/T DNA, without `N` or other ambiguity codes. Keep it in the
  coding orientation. Its length must be divisible by three, and translation with
  the correct genetic code must contain no internal stop codons.
- Do not add or delete arbitrary bases to make the length divisible by three.
  Resolve partial codons or incomplete mitochondrial stop annotations using the
  source annotation and a validated preparation, not automatic padding. A genuine
  terminal stop must agree with the annotated complete CDS; a barcode fragment
  does not need a start or terminal stop codon.

For the usual insect mitochondrial COI analysis, use genetic code **5**. Confirm
the correct code for other organisms. If you cannot verify the interval or its
reading frame, obtain a prepared reference from the maintainer before proceeding.

#### Save and check the guide

From the project root, create its folder:

```bash
mkdir -p data/postprocessing/coi_refs
```

Save the prepared FASTA as `data/postprocessing/coi_refs/COI_reference.fasta`, or
choose a descriptive name and use that path in every command. Renaming alone does
not remove primers, select the right region, or fix a reading frame. Keep the
original downloaded record separately. Beside the guide, record the accession
with version or BOLD process ID, download date, source link, species, coordinates,
strand, expected coding length, genetic code, primer removal and evidence used to
exclude NUMTs. These notes are important when someone else repeats the analysis.

You can run the pipeline's basic checks before a long analysis:

```bash
python -c "import sys; sys.path.insert(0, 'code'); from kcalign_alignment import validate_reference; name, sequence = validate_reference('data/postprocessing/coi_refs/COI_reference.fasta', table=5); print(name); print('Length:', len(sequence), 'nt'); print('Basic checks passed; verify source, exact COI boundaries and absence of NUMTs separately.')"
```

This checks the number of sequences, letters, length divisibility and internal
stops. It does **not** verify species identity, exact ORF boundaries, primer
absence or mitochondrial origin. Compare the reported length with the expected
length established from the source record. Check the alignment after the run as
well; kc-align is not a general guarantee that every primer has been removed.

### A.7 Run COI Postprocessing for one taxon and one barcode

The example uses **Ceratopogonidae** (biting midges) from **barcode01**. Replace
the taxon and barcode if needed. The coding guide must be appropriate to your
group and COI region. Confirm the BOLD database and guide files exist first:

```bash
ls data/bold/BOLD_filtered.fasta.gz
ls data/postprocessing/coi_refs/COI_reference.fasta
```
```bash
python code/postprocess_sintax.py --bold \
  --bold-db data/bold/BOLD_filtered.fasta.gz \
  --taxon Ceratopogonidae --folders barcode01 \
  --analysis-dir data/sintax_analysis/mergedseq5_20261021 \
  --output-dir data/postprocessing/Ceratopogonidae_barcode01_A_run01 \
  --kcalign --kcalign-reference data/postprocessing/coi_refs/COI_reference.fasta \
  --translation-table 5 --threads 6
```

`--folders` selects sequencing barcode folders, not a forward primer such as
`bc1009`. All extracted forward-primer groups in the selected folder are processed.
Put multiword taxon names in quotes, for example `--taxon "Culicoides punctatus"`.
For several barcodes, list them after `--folders`, for example
`--folders barcode01 barcode02`. To use all analysed samples, replace
`--folders barcode01` with the parent argument `mergedseq5`; do not use both.
An `unclassified` folder is allowed only if it was included in the completed analysis.

Use 1 thread or a multiple of three for kc-align, such as **6**. Leave the default
COI protein screen on. Module 2 creates the final results folder; do not create
`data/postprocessing/Ceratopogonidae_barcode01_A_run01` yourself. Choose a new folder name for each new analysis. Copy each whole
command block; a final backslash continues the same command on the next line.

### A.8 Check completion or continue an interrupted run

Keep Ubuntu open. The tools print progress and may be quiet for some time.
A successful run ends with `Completed:` and creates a success marker:

```bash
ls data/postprocessing/Ceratopogonidae_barcode01_A_run01/postprocessing.complete
ls data/postprocessing/Ceratopogonidae_barcode01_A_run01/*.figure.pdf
```

If there is an error, read the last terminal lines or the log:

```bash
tail -n 40 data/postprocessing/Ceratopogonidae_barcode01_A_run01/postprocessing.log
```

Fix the reported problem, make sure no other job is using that results folder,
then repeat the **same** command with `--resume`. For this route:

```bash
python code/postprocess_sintax.py --bold \
  --bold-db data/bold/BOLD_filtered.fasta.gz \
  --taxon Ceratopogonidae --folders barcode01 \
  --analysis-dir data/sintax_analysis/mergedseq5_20261021 \
  --output-dir data/postprocessing/Ceratopogonidae_barcode01_A_run01 \
  --kcalign --kcalign-reference data/postprocessing/coi_refs/COI_reference.fasta \
  --translation-table 5 --threads 6 --resume
```

Keep the same taxon, barcodes, database, guide and results folder for a resume.
Do not create success-marker files yourself. For changed inputs or a new analysis,
use a new results folder rather than resuming the old run.

**Your chosen route ends here. Continue to Section C to understand and review
the results. Do not start the other route.**


For shared file explanations, filtering decisions, UGENE checks and Peartree
viewing, continue to [Section C](post_processing.md#common-results)
or the [HTML results section](metabarcoding-guide.html#review-alignment).
