# Postprocessing (Module 2): choose one starting route

Module 2 makes a tree from already classified DNA reads and named comparison
sequences. **Choose either Route A or Route B. You do not need to follow both.**
Both routes below focus on one taxon, one barcode, BOLD comparison sequences and
kc-align for COI. MAFFT remains the software default, but kc-align is the preferred
COI choice when you have a trustworthy coding guide.

| Your starting point | Follow | Do not repeat |
| --- | --- | --- |
| You downloaded `postprocessing-inputs.tar.gz` | **Route A: shared inputs** | Module 1 and Route B. The archive already contains the processed inputs. |
| You processed your own reads in Module 1, or need to do so | **Route B: your own analysis** | Archive extraction, archive path relocation and Route A. |
| Either route has finished | **Section C: understand and review your results** | No second Postprocessing run is required just to read the results. |

Each route is complete: database and guide preparation are repeated deliberately
so you can stay in your chosen route. Computer installation is covered in
[setup](setup.md). The [HTML guide](metabarcoding-guide.html#postprocessing)
provides the same routes with colour-coded navigation and copy buttons.

<a id="if-you-downloaded-postprocessing-inputstargz-start-here"></a>
## Route A — downloaded postprocessing-inputs.tar.gz

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

## Route B — your own Module 1 analysis, without the input archive

**Use only this route for your own sequencing data.** You will use the completed
Module 1 analysis produced on your computer, or create it if it is not ready.
**Do not download or extract the shared archive. Do not run Route A's relocation
block on your own analysis.** Once Module 1 is complete, you do not run it again
for every taxon or barcode.

### B.1 Choose a stable Linux project location

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

If Module 1 is already complete in a Windows folder, the simplest way to keep
its saved paths valid is to continue using that project at its existing location:

```bash
cd "/mnt/c/Users/USERNAME/Documents/metabarcoding"
pwd
```

This is supported, though it may be slower. If you want to move a completed
analysis to Linux, obtain a path-relocation procedure for **your own dataset**
from the maintainer; Route A's example targets a different saved analysis.
For a new analysis, choose Linux before starting Module 1 and keep that location.

### B.2 Activate the software

From the project root, create the environment **only if it has not been
installed**:

```bash
conda env create --name taxonomy --file config/taxonomy-environment.yml
conda activate taxonomy
```

If it already exists, run only `conda activate taxonomy`. Activate it whenever you
open a new terminal for analysis. If Git or Conda is not available, finish the
[computer setup instructions](setup.md) before continuing.

### B.3 Finish Module 1 only if it is not already complete

If you already have a completed analysis, skip the processing commands in this
step and go directly to B.4. Otherwise put your barcode folders under your own
dataset parent, for example `data/sequencing2/`. Replace `sequencing2` consistently
if yours has a different name. See [Read processing](pipeline_from_fastplong_onwards.md)
for required laboratory primer/adapter files and classifier-database setup.

For a new Linux analysis, a dataset still in Windows can be copied before
processing, replacing the example source path:

```bash
cd ~/metabarcoding
mkdir -p data
cp -a "/mnt/c/Users/USERNAME/Documents/sequencing2" data/sequencing2
```

Run that copy only if `data/sequencing2` does not already exist, to avoid nesting
or mixing datasets. It leaves the Windows original unchanged. `cd` alone would
not copy those files.

If each barcode contains nested `pass/` and `fail/` FASTQs, merge them first:

```bash
python code/merge_pass_fail.py sequencing2 --dry-run
python code/merge_pass_fail.py sequencing2
```

This keeps the original files and creates a combined gzip file per sample.
Skip the merge if your files are already directly inside each barcode folder.
Then run Module 1 once:

```bash
python code/read_processing.py sequencing2 --output-dir data/read_processing/sequencing2_run01
```

Wait for Module 1 to finish. Do not start Module 2 from a folder containing only
raw reads. Do not rerun Module 1 just because you want another taxon later.

### B.4 Identify your completed analysis and available samples

If you used B.3's output path, the completed analysis is
`data/read_processing/sequencing2_run01/analysis`. Use your actual analysis path
if it differs; do not choose its surrounding run folder or `primertrack/`.

```bash
ls data/read_processing/sequencing2_run01/analysis/analysis.complete
python -c "import json; p=json.load(open('data/read_processing/sequencing2_run01/analysis/analysis_parameters.json')); print('Dataset:', p['merged_dir']); print('Samples:', ', '.join(p['analysed_folders']))"
```

The first check must find the marker. The second lists the dataset and sample
names. Select samples from that list. Keep the original cleaned reads, saved
classification TSVs and completed analysis at their recorded locations. Module 2
uses those results; it does not rerun the classifier or open its source database.

### B.5 Download and prepare BOLD

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


Use `data/bold/BOLD_filtered.fasta.gz` in B.7. Choosing that one file explicitly
avoids accidentally reading old or duplicate snapshots. The preparation removes
empty records, while Module 2 later retains named-species references and applies
its 500–2,000 nucleotide length filter.

### B.6 Choose and prepare the COI coding guide for kc-align

The guide is a separate, single-sequence file. The database and the Module 1
analysis do not automatically prepare it. These full selection and checking
instructions are repeated here so you can stay in Route B.

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

### B.7 Run COI Postprocessing for one taxon and one barcode

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
  --analysis-dir data/read_processing/sequencing2_run01/analysis \
  --output-dir data/postprocessing/Ceratopogonidae_barcode01_B_run01 \
  --kcalign --kcalign-reference data/postprocessing/coi_refs/COI_reference.fasta \
  --translation-table 5 --threads 6
```

`--folders` selects sequencing barcode folders, not a forward primer such as
`bc1009`. All extracted forward-primer groups in the selected folder are processed.
Put multiword taxon names in quotes, for example `--taxon "Culicoides punctatus"`.
For several barcodes, list them after `--folders`, for example
`--folders barcode01 barcode02`. To use all analysed samples, replace
`--folders barcode01` with the parent argument `sequencing2`; do not use both.
An `unclassified` folder is allowed only if it was included in the completed analysis.

Use 1 thread or a multiple of three for kc-align, such as **6**. Leave the default
COI protein screen on. Module 2 creates the final results folder; do not create
`data/postprocessing/Ceratopogonidae_barcode01_B_run01` yourself. Choose a new folder name for each new analysis. Copy each whole
command block; a final backslash continues the same command on the next line.

### B.8 Check completion or continue an interrupted run

Keep Ubuntu open. The tools print progress and may be quiet for some time.
A successful run ends with `Completed:` and creates a success marker:

```bash
ls data/postprocessing/Ceratopogonidae_barcode01_B_run01/postprocessing.complete
ls data/postprocessing/Ceratopogonidae_barcode01_B_run01/*.figure.pdf
```

If there is an error, read the last terminal lines or the log:

```bash
tail -n 40 data/postprocessing/Ceratopogonidae_barcode01_B_run01/postprocessing.log
```

Fix the reported problem, make sure no other job is using that results folder,
then repeat the **same** command with `--resume`. For this route:

```bash
python code/postprocess_sintax.py --bold \
  --bold-db data/bold/BOLD_filtered.fasta.gz \
  --taxon Ceratopogonidae --folders barcode01 \
  --analysis-dir data/read_processing/sequencing2_run01/analysis \
  --output-dir data/postprocessing/Ceratopogonidae_barcode01_B_run01 \
  --kcalign --kcalign-reference data/postprocessing/coi_refs/COI_reference.fasta \
  --translation-table 5 --threads 6 --resume
```

Keep the same taxon, barcodes, database, guide and results folder for a resume.
Do not create success-marker files yourself. For changed inputs or a new analysis,
use a new results folder rather than resuming the old run.

**Your chosen route ends here. Continue to Section C to understand and review
the results. Do not start the other route.**

<a id="common-results"></a>
## Section C — understand, inspect and decide what to keep

**This section is common to both routes.** It explains files and quality checks;
it is not another set of extraction or Postprocessing steps. Do not rerun the
other route. Your analysis folder is the path you supplied after `--output-dir`.

### C.1 Find the important files

For Route A's example, open `data/postprocessing/Ceratopogonidae_barcode01_A_run01/`.
For Route B's example, open `data/postprocessing/Ceratopogonidae_barcode01_B_run01/`.
The `*` below stands for the varying taxon/date prefix, not a literal filename.

| File | Produced by | Contents and purpose |
| --- | --- | --- |
| `postprocessing.complete` | Module 2 | Success marker; a picture alone does not confirm the latest run finished. |
| `*.figure.pdf`, `*.figure.png`, `*.figure.svg` | Baltic plotting helper | Tree pictures for reading and sharing. |
| `*.figure.midpoint.nwk` | Baltic plotting helper | Rooted Newick tree to open in Peartree. |
| `extraction_by_barcode_primer_taxon.tsv` | Module 2 from local extraction | Counts, forward/reverse primers, taxonomic labels and retained/excluded decisions. |
| `cleaned_*.fasta` | Extraction helper | Original extracted local reads; `.dominant.fasta` keeps selected reads and `.vsearch97.fasta` contains cluster consensuses. |
| `vsearch_local_groups.tsv`, `vsearch_local_memberships.tsv` | Module 2 / VSEARCH | Local group counts and the reads contributing to each cluster. |
| `vsearch_reference_groups.tsv`, `vsearch_reference_memberships.tsv` | Module 2 / VSEARCH | Reference taxon groups and cluster membership. |
| `bold_reference_metadata.tsv` | BOLD extraction helper | Original BOLD IDs, labels and source headers for selected references. |
| `*.aln.fasta` | kc-align or MAFFT | Untrimmed DNA alignment; keep for comparison. |
| `*.clipkit.aln.fasta` | ClipKIT | Trimmed alignment used for tree inference; inspect this in UGENE. |
| `*.kcalign_selection.tsv` | kc-align adapter | Which consensus sequences were aligned or rejected. |
| `*.treefile`, `*.iqtree` | IQ-TREE | Estimated tree and model/inference report; new runs have no bootstrap replicates. |
| `*.figure.tips.tsv` | Plotting helper | Links tree labels to local reads or database records. |
| `parameters.json`, `checkpoints.json`, `postprocessing.log` | Module 2 | Settings, stage/output checks and commands/errors for auditing or retrying. |

Open `.tsv` files in a spreadsheet and choose **Tab** as the separator if asked.
For descriptions of all main outputs, see [the output-file reference](../data/postprocessing/README.md).

#### Read the local FASTA headers: a link back to your data

A FASTA header is the line beginning with `>` above each DNA sequence. Module 2
adds tracking information to that line. The current local-read layout is:

```text
>barcode-forward_primer-reverse_primer-trusted_taxon-contamination_flag-read_id
```

For example, this invented header illustrates the structure:

```text
>barcode01-bc1009-none-Culicoides_punctatus_12345-nocontamination-read-abc-001
```

| Part | Meaning | Where to follow it |
| --- | --- | --- |
| `barcode01` | Source sequencing sample folder. | The corresponding Module 1 sample read table and cleaned `filtered/asinput.fastq.gz`. |
| `bc1009` | Canonical forward-primer annotation, obtained from tracking the raw parent read. | Module 1's `forward_PCR_primer` field and Module 2's extraction-count table. |
| `none` | Reverse-primer annotation; here a tracking record exists but no reverse primer was detected. | Module 1's `reverse_PCR_primer` and Module 2's `reverse_primer` column. `unknown` instead means no matching tracking record. |
| `Culicoides_punctatus_12345` | Lowest taxonomic label trusted at the saved confidence cutoff, normally 0.80. The final numeric part, when present, is a taxonomy ID—not a read count. | Module 1's `assigned_label`/`assigned_name` and saved classification TSV, then `taxonomic_identification` in the extraction table. `12345` here is only an example. |
| `nocontamination` | The saved analysis did not flag that trusted label by its across-barcode prevalence rule. It is not a guarantee that this particular read is contamination-free. | The saved taxon-prevalence report. `contamination` is the alternative flag; this flag is distinct from the highest-read-count selection. |
| `read-abc-001` | Original cleaned-query identifier, which can itself contain hyphens. | Module 1's `readid` column and the sequence ID in cleaned FASTQ/SINTAX outputs. `raw_ont_read_id` in the read table links it to the original parent when adapters split a read. |

The **first five hyphens** separate the six fields; the rest of the text is the
read ID. Header annotations escape hyphens inside their own values, for example
`%2D`, so those are not mistaken for separators. Spaces in annotations become
underscores, and other special characters may appear as percent escapes. Do not
rename or split headers arbitrarily: extraction, clustering, alignments and tip
labels use these identifiers to keep the records connected.

Follow a sequence through the outputs in this order:

1. In `cleaned_TAXON_DATE.fasta`, the header describes an extracted original local
   read. The extraction TSV counts reads by barcode, forward primer, reverse
   primer and trusted taxon. Winning-label selection sums across reverse primers
   for each barcode + forward primer; all tied highest-count labels are kept.
2. In `.dominant.fasta`, retained reads keep their original tracking headers.
   The prevalence flag does not decide which label wins this count comparison.
3. In `.vsearch97.fasta`, the sequence is a **cluster consensus**, but its header
   is retained from the seed read. Use `vsearch_local_groups.tsv` and
   `vsearch_local_memberships.tsv` to find the group and all contributing reads.
   Reverse primers and contamination flags in that seed's name need not describe
   every member; they do not create additional local cluster groups.
4. The header carries through the alignment and ClipKIT-trimmed alignment to the
   tree. Translation filtering or kc-align can remove a consensus; check
   `translation_lengths.tsv` or `*.kcalign_selection.tsv` when a header disappears.
5. `*.figure.tips.tsv` connects a final tree tip to the seed read, source folder
   and display label. A tip name is tracking information, not a fresh taxonomic
   classification of the consensus and not necessarily a separate insect.

BOLD and NCBI references have a different layout, for example
`BOLD_000000123__taxon__Culicoides_punctatus` or
`AB123456.1__taxon__Culicoides_punctatus`. The part after `__taxon__` is the safe
species/organism label used to keep reference taxa separate during clustering.
For BOLD, use `bold_reference_metadata.tsv` to recover the original process ID,
full header and database source. For NCBI, use `reference_records.gb.xml` and the
accession. These references do not have local barcode or primer annotations.


### C.2 Understand which reads are kept



The program performs these steps automatically:

1. **Extract your reads.** Select the requested organism group from the chosen
   barcode folders, using the identifications recorded by Module 1.
2. **Keep the highest-count identification per barcode + forward primer.** The
   assumption is one insect for each such pair. Counts are added across reverse
   primers. For example, counts of 37 for Tipula tipula, 10 for Tipula spectra and
   2 for Tipula nigris retain the 37 Tipula tipula reads. If labels tie for the
   largest total, all tied labels and their reads are retained. Lower-count
   labels do not continue. This is calculated only among the extracted reads;
   the winner need not exceed half the total and is not proof of species identity.
3. **Group similar sequences.** Local reads cluster at **97% identity**, separately
   by barcode, forward primer and retained label. Each cluster produces a
   *consensus*: a summary DNA sequence built from its reads. A broad genus-level
   assignment can hide different species, so sequence similarity is still checked.
   Reference sequences cluster at **98%**, only within the same taxonomic label;
   different reference species never share a cluster. Local and reference reads
   are kept separate during clustering.
4. **Check protein-coding DNA.** For the usual insect COI analysis, keep the default
   protein check. It requires at least 200 amino acids in a row without a stop or
   unknown letter in one possible reading frame, using genetic code 5. Passing
   does not prove that the sequence is genuine mitochondrial COI.
5. **Align and trim.** For COI, prefer kc-align with a prepared coding guide;
   MAFFT runs if you omit that option. The aligner lines up corresponding DNA. ClipKIT
   removes columns with many gaps. Both steps are automatic.
6. **Estimate and draw the tree.** IQ-TREE estimates the tree and Baltic draws it.
   Bootstrap replicates are disabled, so new trees do not have bootstrap support.
   The picture is midpoint rooted: its root is placed halfway along the longest
   path between tips. This does not establish the organisms' biological ancestor.

**Primer removal needs attention:** tracking a primer is not the same as removing
it. The current MAFFT workflow has no general primer-removal step. If primer
sequences remain in your reads, arrange a suitable removal workflow or have the
alignment reviewed before interpreting the resulting tree. Removing gap-heavy
columns is not a substitute for primer removal.

**Noncoding DNA:** use MAFFT with `--skip-translation` and references for the same
marker. Do not change the genetic code or protein threshold just to make a failed
run finish. Confirm the marker first.

**Preferred for COI: kc-align.** use only with a prepared, trusted coding reference for your
marker. Add `--kcalign --kcalign-reference PATH`. The FASTA must contain exactly one
ungapped A/C/G/T sequence in the correct reading frame, without internal stops;
its length must be divisible by three and it should not contain primer sequence.
If you cannot verify this, obtain a prepared reference from the maintainer.
Use 1 thread or a multiple of three, such as 6. It can reject reads; its selection
TSV records those decisions. It does not automatically add the guide as a tree tip.

<a id="review-the-trimmed-alignment-before-interpreting-the-tree"></a>

### C.3 Check the trimmed alignment in UGENE and choose records



**You must check the alignment after ClipKIT.** Finishing the pipeline is not a
substitute for reviewing the DNA comparison. The current pipeline proceeds from
ClipKIT to IQ-TREE automatically; it does not pause for manual approval. Treat
its first tree as a draft until you have checked the alignment and inputs.

#### Open the alignment in UGENE

[Download UGENE](https://ugene.net/download-all.html), a desktop program for
viewing and editing sequence alignments. Choose the installer or package for
your operating system. On Windows, you can install and run UGENE normally while
the pipeline runs in Ubuntu; it does not need to be installed in the Conda
environment. See the [UGENE alignment-editor guide](https://ugene.net/quick-start-guide/what-you-can-do-in-ugene/view-edit-and-align-multiple-sequence-alignments/).

1. Find your run's **`*.clipkit.aln.fasta`** file: this is the trimmed alignment
   used for its tree. Also keep the corresponding untrimmed `*.aln.fasta` file
   available for comparison.
2. In UGENE, choose **File → Open** and select the trimmed file. If asked how to
   open a multi-sequence FASTA, choose to join/open it as an alignment, preserving
   the existing gaps; do not join the sequences end to end or realign them merely
   to view the file.
3. For a Linux project viewed from Windows, browse through `\\wsl$`, your Ubuntu
   distribution, `home`, your Linux username and `metabarcoding/data/postprocessing`.
   A project under `/mnt/c/` is accessible through its usual Windows folder.
4. Scroll across the entire alignment and compare local consensuses with trusted
   references. Check sequence ends, possible remaining primers, blocks that do
   not line up, many gaps or ambiguous letters, unusually short sequences,
   unexpected insertions/deletions and suspiciously divergent records. For COI,
   assess the coding alignment against the trusted guide and its translation;
   ClipKIT's nucleotide column trimming does not itself establish a correct ORF.
5. Decide whether the sequence set and aligned region are suitable. Keep the
   original outputs and save any edited alignment as a separate reviewed copy
   outside the pipeline's output folder, with notes explaining every change.

If you remove sequences or change alignment columns, **estimate the tree again
from the reviewed alignment** before interpreting it. Hiding a tip in a viewer
does not re-estimate the relationships. Ordinary `--resume` is not a manual
alignment-import feature: it validates saved outputs and may rebuild changed
files. Do not overwrite its alignment or success markers and expect a new tree
from your edits. Use a documented separate tree inference from the reviewed
alignment, or prepare a revised input dataset and start a new pipeline analysis.

#### Keep numbered alignment and tree versions

In UGENE, save/export each edited alignment under a **new filename**, rather than
saving over the ClipKIT output or your previous reviewed copy. Use a simple
version number and keep the original sequence headers unchanged. For example:

| File or IQ-TREE prefix | Meaning |
| --- | --- |
| `data/final/coi_aln_v01.fasta` | First reviewed alignment. |
| `data/final/coi_aln_v02.fasta` | Next alignment version after further edits. Keep v01 too. |
| `data/final/coi_aln_v01_run01` | IQ-TREE output prefix for the first estimation using alignment v01. |
| `data/final/coi_aln_v02_run01` | A separate estimation using alignment v02. |
| `data/final/coi_aln_v02_run02` | Another estimation from v02, for example with a different model or seed. |

Add the taxon/barcode or date to the names if you are managing several datasets.
Keep a short `data/final/README.md` recording each alignment version, what changed,
which reads or reference species were removed and why, and the exact IQ-TREE
command used. This links every tree to the alignment that produced it.

**Use a new `--prefix` for each new phylogenetic estimation you want to keep.**
IQ-TREE writes several files with that prefix, including `.treefile`, `.iqtree`
and `.log`; changing only the displayed picture name does not preserve a separate
inference. Do not use `--redo` to overwrite a run you want to retain. Continuing
an interrupted IQ-TREE job with the same input and command is different from
starting a new estimation and can reuse that job's prefix.


#### Estimate the final reviewed phylogeny with bootstrap support

The automatic Module 2 run disables bootstrap to make the initial analysis faster.
**Once you have finished evaluating and editing the alignment and chosen the
sequences to retain, use bootstrap for your final IQ-TREE phylogeny.** This is a
separate final inference, not a change to Module 2's quick-run default.

Create a separate final-results folder and save/export your reviewed alignment
there as `coi_aln_v01.fasta`, preserving its FASTA headers and existing aligned
columns. If you made no edits, use a copy of the checked ClipKIT alignment:

```bash
mkdir -p data/final
```

After placing the reviewed file at `data/final/coi_aln_v01.fasta`, run:

```bash
iqtree -s data/final/coi_aln_v01.fasta -m MFP \
  -B 1000 --bnni -T 6 --seed 666 \
  --prefix data/final/coi_aln_v01_run01
```

This example uses **1,000 ultrafast bootstrap replicates** (`-B 1000`) with
additional bootstrap-tree optimisation (`--bnni`), not traditional nonparametric
bootstrap. It takes longer than the initial run. Use a new output prefix for a
different final analysis. Check `coi_aln_v01_run01.iqtree` for the model, convergence
information and support method; inspect `coi_aln_v01_run01.treefile` or the bootstrap
consensus `coi_aln_v01_run01.contree`. Support values help assess branch stability but do
not establish correct species identity or compensate for a poor alignment.
See the [IQ-TREE bootstrap documentation](https://iqtree.github.io/doc/Command-Reference#ultrafast-bootstrap-parameters).

These final outputs are separate from the pipeline's original pictures and
checkpoint files. Load the new tree when viewing the final result; the earlier
PDF does not automatically update. Ordinary Module 2 `--resume` will not import
your manually reviewed alignment or add this bootstrap analysis.


#### Which reads and reference species should stay?

**The final choices are at the user's discretion and scientific judgement.**
You decide which local reads or consensuses to exclude, which BOLD species and
individual reference records to retain, and what comparison set is appropriate
for the question. The pipeline's highest-count-label, species, length and
clustering filters provide a starting set; they do not make those final
biological decisions for you.

Base exclusions on evidence such as poor sequence quality, primer contamination,
questionable identifications, suspected NUMTs, inadequate overlap or an unsuitable
marker. Review the extraction table, cluster memberships, database record and
alignment together. Record IDs, source releases and reasons so another person
can repeat the selection. Keep the unmodified inputs. Do not just change
`kept_for_module2` in the generated table: it is an audit of automatic selection,
not a supported control for manually choosing pipeline inputs.

### C.4 Explore the tree in Peartree



Open [Peartree](https://peartree.live/) in your browser and load the run's
**`*.figure.midpoint.nwk`** file. This is a Newick tree with the midpoint root used
in the pipeline's figure. Load the tree file, not the PDF or PNG image. You can
also inspect the IQ-TREE `*.treefile`, but it may be displayed with a different
root. See the [Peartree manual](https://peartree.live/manual/) for viewing controls.

Use Peartree to explore the tree and read crowded labels. A viewer's display,
rooting or hidden-tip changes do not change the alignment or rerun inference.
The PDF, PNG and SVG supplied by the pipeline remain available for quick viewing
and sharing. Check the alignment in UGENE as well; viewing the tree alone is not
an alignment-quality check. Bootstrap support remains disabled in new runs.



A tip is a cluster consensus, not necessarily one raw read or one insect. Its
header describes the cluster's seed, not a fresh classification of the consensus.
Branch lengths represent DNA change, not years. Midpoint rooting is a plotting
choice, not proof of the biological ancestor. Bootstrap replicates are disabled.
A local tip near a named reference is evidence for comparison, not a confirmed
species identification on its own. Your read, reference-species and alignment
choices remain your scientific decisions; document them and rebuild the tree
when its input alignment changes.

#### Optional Python figures with Baltic

For custom tree visualisation in Python, **Baltic** is recommended. The pipeline
already uses it for its tree figures, and it is included in the `taxonomy`
environment. Browse the [Baltic gallery, tutorials and examples](https://phylo-baltic.github.io/baltic-gallery/)
for ways to colour tips, add metadata and customise figures. This is optional;
Peartree remains a quick way to explore the tree without writing Python code.

For a final bootstrap tree, load the new `data/final/coi_aln_v01_run01.treefile` or
`.contree` in Baltic or Peartree, rather than the initial pipeline tree. Preserve
and label the support method when plotting. The project's
[plotting helper](../code/plot_phylogeny.py) can also render an existing tree with
its corresponding analysis and reference metadata; use `--help` to see the inputs.
Changing the picture's rooting or layout does not repeat the inference.


### C.5 Alternatives and help

For noncoding markers, use MAFFT with `--skip-translation` and matching references;
do not apply COI coding-guide advice to a noncoding marker. For NCBI comparison
references, replace `--bold` with `--download` and omit `--bold-db ...`; NCBI mode
searches COI and needs internet and Entrez Direct. Its download does not supply
the kc-align guide. MAFFT remains the command-line default if `--kcalign` is absent.

If something fails, first read the error and the last lines of your run's log.
Check input paths, environment activation, taxon spelling, sample availability,
reference quality and whether at least four aligned sequences remain. Send the
full command and error to the maintainer rather than changing validation values.
Use the optional command builder in the HTML if you want a different sample set;
it does not run the analysis or check that files exist.

## Appendix — optional technical reference

The routes above are complete. The following sections document tool internals;
they are not a third workflow. Their historical example paths are not defaults.


The walkthrough above is enough for a standard run. The remaining sections
retain the exact commands, file fields and specialist settings for readers who
need them. Example paths below may refer to older project datasets; substitute
your own completed analysis and input paths.


For Linux Miniconda installation, accessing Windows folders from Ubuntu, and
cloning the [GitHub repository](https://github.com/mftorres/metabarcoding), follow
[Ubuntu/WSL setup](setup.md).

[Read processing](pipeline_from_fastplong_onwards.md) cleans reads, assigns taxonomy with
SINTAX and generates read-level analyses. Postprocessing exports selected cleaned reads,
obtains external references, builds cluster consensus sequences, optionally screens
translations, aligns nucleotides, estimates an IQ-TREE phylogeny and draws a Baltic
midpoint-rooted divergence tree. Run from the repository root in `taxonomy`.

For a size-conscious collaborator upload, see [the required-file lists and tar.gz instructions](sharing_postprocessing.md).

### Run the module

```bash
conda activate taxonomy
python code/postprocess_sintax.py sequencing2 --bold \
  --taxon Tipulidae \
  --analysis-dir data/read_processing/sequencing2_run01/analysis \
  --output-dir data/postprocessing/Tipulidae_bold_run02 --threads 8
```

Pass the parent once (`sequencing2`, `data/sequencing2`, or a full path): the script
selects its folders from the completed Read processing analysis, including `unclassified`
when that folder was analysed. `--folders barcode01 barcode02` remains available
for a subset, and `--folders data/sequencing2` also accepts a parent. Omit both to
select all analysed folders. A supplied parent must match the dataset saved in
the analysis metadata. Do not combine a positional parent with `--folders`.

For the [complete parent-folder workflow and pass/fail merge commands](pipeline_from_fastplong_onwards.md#prepare-passfail-folders-before-analysis),
run `python code/read_processing.py sequencing2 --output-dir data/read_processing/sequencing2_run01`
first. The launcher prepares the `analysis/` directory used above. Unclassified
reads have their own source folder and are excluded from barcode prevalence in Read processing;
they remain traceable as local tips in Postprocessing. Other taxon, source and aligner
options are still required as appropriate.

Use a new output directory for a new analysis. Dependencies are `vsearch`, `seqkit`,
`clipkit`, `iqtree`, Python with Baltic >=1.0.0 and Matplotlib, and the chosen aligner: `mafft`
(default) or `kc-align`. NCBI mode additionally requires Entrez Direct `esearch` and
`efetch`. kc-align is installed in the `taxonomy` environment.

### Extracted-read abundance table and output reference

Local extraction first keeps all requested reads in `cleaned_TAXON_DATE.fasta`.
The audit table `extraction_by_barcode_primer_taxon.tsv` has columns `barcode`,
`forward_primer`, `reverse_primer`, `taxonomic_identification`, `read_count` and
`kept_for_module2` (`true`/`false`). Counts are split by reverse primer, but selection
sums each taxonomic label across reverse primers within each barcode + forward
primer. Every label tied for the largest total survives. All reads with any of those
labels are retained regardless
of reverse primer or contamination flag; other labels are excluded from downstream
analyses. Selection uses only the requested taxon's extracted reads, not reads
outside the extraction. `none`, `unknown` and combined primer annotations remain
exact distinct forward-primer values, including in `unclassified` when requested.

The assumption is one insect individual per barcode + forward primer, whose DNA
should dominate. The most abundant label need not exceed 50% and is not proof of
identity or absence of contamination. For example, barcode15/bc1009 with Tipula
tipula (37), Tipula spectra (10) and Tipula nigris (2) retains the 37 Tipula tipula
reads. `cleaned_TAXON_DATE.dominant.fasta` contains the retained original reads;
the original extraction and every audit-table row remain available.

Counts represent cleaned queries, including separate adapter-split fragments,
not distinct raw parents or final tree tips. Reverse-primer annotations in a cluster
header describe its seed; contributing reads may have different reverse primers.
A compatible resume regenerates the table when the original extraction exists.

For the contents, producing program and interpretation of every main output,
see [the Postprocessing folder README](../data/postprocessing/README.md), including
the file patterns found in `Ceratopogonidae_bold_kcalign_run01/`.

### Choose the reference dataset

For local BOLD references, follow the [database download and preparation guide](../data/bold/README.md).
It links to BOLD's Data Packages page, explains where to put the gzip FASTA
(`data/bold/`), and shows `seqkit seq --min-len 1` to remove records without
sequences. BOLD stays compressed; Postprocessing reads gzip directly. Keep original
unfiltered downloads outside the scanned database folder.

Exactly one source flag is required:

- `--download`: resolve a name or numeric taxid with NCBI, then download COI/CO1/COX1
  nucleotide records including descendants. Searches also include cytochrome
  oxidase I/1 and cytochrome c oxidase subunit I/1 title aliases.
- `--bold`: stream `.fasta.gz`, `.fa.gz` and `.fna.gz` databases in `data/bold/`
  locally. `--bold-db PATH` selects one file or another directory. This mode accepts
  a taxon name, not a numeric taxid, and makes no Entrez calls. It matches a substring
  of the full FASTA header, ignoring case and treating underscores as spaces.
  Every matching marker is included; it does not implicitly restrict markers to COI.

BOLD headers are expected to contain `processID|marker|country|comma-separated lineage`.
The original header stays in the FASTA description. Safe IDs `BOLD_000000001`, etc.,
are mapped to process IDs and taxonomy by `bold_reference_metadata.tsv` with columns
`sequence_id`, `process_id`, `marker`, `taxon`, `original_header`, `source_file`.
The display taxon is the species identification; genus-only and provisional records are excluded.

Mitogenomes are excluded by title by default; `--include-mitogenomes` permits them.
This title-based rule is not a guarantee that every mitogenome is detected.
Both sources are then filtered with `seqkit seq --min-len 500 --max-len 2000`:
the bounds are inclusive. This reference length filter remains active even when
translation is skipped. Empty reference datasets stop the run.

### Local read selection and headers

`export_sintax_taxa.py` selects reads when the queried name, complete classifier label
or taxid matches a rank in the full SINTAX lineage **and that rank has confidence
>=0.80**. Lower ranks need not pass. Tipulidae at 0.90 with Tipula at 0.20 qualifies
for a Tipulidae query but not a Tipula query. Missing confidence does not qualify.
SINTAX confidence is bootstrap support, not a calibrated probability.

Cleaned sequences are extracted by read ID. Raw reads are not used as phylogenetic
input. Headers contain:

```text
folder-forwardprimer-reverseprimer-lowest_trusted_taxon-contamination_flag-readID
```

The taxon is the lowest rank with confidence >=0.80; contamination applies to that
trusted taxon according to the completed analysis. Primer names use the canonical
mapping in that analysis. Header components escape delimiter characters; read IDs
retain their tracking information. Extraction counts for every requested folder
and the total are saved in `extraction_counts.tsv` and `extraction.log`.

### Clustering: 97% local and 98% reference consensus sequences

After dominant-taxon selection, local reads are clustered independently for each
**barcode + forward primer + retained trusted taxonomic label**. Reverse primers
are tracked in the audit but do not split an individual into separate clusters.
VSEARCH still applies a **97% identity cutoff** within each individual. A retained
genus or higher-rank label can hide multiple species without references; the
identity cutoff helps keep their reads in separate consensuses rather than
combining all reads with the same assignment into a single consensus.

`none` means a primer-tracking record exists for the read's raw parent, but no
primer was detected in that direction. `unknown` means no matching primer-tracking
record was found, so the primer assignment is unavailable. A missing record is
not evidence that the primer was absent. These meanings apply separately to the
forward and reverse primer fields.

BOLD extraction retains only species-identified records. Both BOLD and NCBI references are partitioned by
the exact taxonomic suffix after `__taxon__` in their FASTA headers, then clustered
at 98% within each label. Different species labels never share a consensus, even
if their sequences are identical. This preserves species positions in the tree
instead of hiding them under a different species' seed label. NCBI labels remain exact groups; no new species identification is inferred. References
and local reads remain separate. Local groups use the command below; reference
groups use the same command with `--id 0.98`:

```bash
vsearch --cluster_fast input.fasta --id 0.97 --iddef 2 \
  --strand both --qmask none --notrunclabels \
  --consout consensus.fasta --uc clusters.uc --threads 8
```

The module sets minimum length 1 and maximum length to the longest input so this
step adds no length filter. Identity definition 2 excludes terminal gaps.
Downstream sequences are **consensuses**, potentially different from individual
input reads. VSEARCH's `centroid=` and `;seqs=N` label decorations are removed;
the seed header is retained for tracking. Its taxon, contamination and read ID or
accession describe the **seed**, not a new classification of the consensus.
No new SINTAX classification is performed. UC files identify all cluster members.

Outputs include:

- `cleaned_TAXON_DATE.vsearch97.fasta`: combined local consensus sequences.
- `COI_TAXON_DATE.vsearch98.fasta` or `BOLD_TAXON_DATE.vsearch98.fasta`: reference consensuses.
- `vsearch_local_groups.tsv`: `group_id`, `source_folder`, `forward_primer`,
  `taxon_label`, `input_sequences`, `consensus_sequences`, `identity`.
- `vsearch_local_groups/`: each group's inputs, consensus FASTA, UC and stdout log.
- `vsearch_local_memberships.tsv`: group ID followed by the ten UC fields.
  Cluster numbers are relative to their group.
- `vsearch_reference_groups.tsv`: reference taxon groups and counts.
- `vsearch_reference_memberships.tsv`: reference group ID followed by ten UC fields.
- `vsearch_reference_groups/`: per-taxon input, consensus FASTA, standard `.uc` and log files.
- `cleaned_TAXON_DATE.dominant.fasta`: local reads retained before clustering.
- `clustering_counts.tsv`: `dataset`, `input_sequences`, `consensus_sequences`, `identity`.

The local and reference consensus sets are concatenated to `combined_TAXON_DATE.fasta`.
There is no further clustering after combining them.

### Translation selection and nuclear sequences

By default, SeqKit translates the combined consensuses in all six frames with
`--translation-table 5` (invertebrate mitochondrial). At least one frame must
contain **at least 200 consecutive amino acids**, splitting at stops `*` and
unknown residues `X`. `--min-aa` adjusts the minimum. Outputs are
`combined_TAXON_DATE.six_frames.faa`, `translation_lengths.tsv`, and the retained
nucleotides `combined_TAXON_nonumt_DATE.fasta`. This is a screening heuristic,
not proof of the absence of NUMTs. The original consensus nucleotides are retained.

**`--skip-translation`** bypasses both translation and amino-acid selection and
keeps all consensuses in `combined_TAXON_selected_DATE.fasta`. It creates no
six-frame FASTA or translation-length table. Use this with default MAFFT for
**nuclear loci, particularly noncoding sequences**, for which the mitochondrial
translation screen is inappropriate:

```bash
python code/postprocess_sintax.py --bold --bold-db data/bold/nuclear_marker.fasta.gz \
  --taxon Tipulidae --folders barcode{01..21} \
  --analysis-dir data/sintax_analysis/mergedseq5_20261021 \
  --output-dir data/postprocessing/Tipulidae_nuclear_run01 \
  --skip-translation --threads 8
```

The nuclear database path above is a placeholder for your marker-specific database,
using the BOLD header format. Skipping translation does not change marker selection:
**NCBI `--download` still searches COI**, and BOLD mode includes all matching markers.
The local reads, their SINTAX reference assignments, and the external database must
be appropriate to the locus being analysed. The 500–2,000 nt reference filter remains.

### Alignment choice

MAFFT is the default (`--mafft` may also be explicit). Current arguments are:

```bash
mafft --retree 2 --maxiterate 0 --adjustdirectionaccurately \
  --op 20 --ep 18 --thread 8 selected_consensus.fasta
```

This is FFT-NS-2 nucleotide alignment, suitable for nuclear sequences as well as
other nucleotide loci. Strand-adjusted headers are restored and recorded in
`alignment_orientation.tsv`. All representative IDs, equal alignment lengths and
preservation of nucleotide sequences (allowing reverse complementation) are checked.

**Important:** For now, the pipeline lacks a step for primer sequence removal. Before tree estimation, the user **must** remove the primer sequences from the alignment. In a next version of the pipeline, I will add the primer sequences as "adapter/barcode" to `fastplong` to remove them at early stages. I doubt (but will test) taxonomic classification will be severely affected by the presence of primer sequences in the clean reads (2026-Oct-29).

Add **`--kcalign --kcalign-reference PATH`** for reference-guided mixed-mode
codon alignment. The reference file is required: **provide a nucleotide sequence
you trust for this locus**, in the correct orientation and reading frame. It must
contain exactly one ungapped A/C/G/T sequence. Its length must be divisible by
three, and it must start at the first nucleotide of the first codon in the
reference coding region. Remove incomplete terminal codons using the known
reading frame; do not trim arbitrarily just to make the length divisible by three.
The reference must have no internal stops under the chosen genetic code. Automatic validation cannot establish
that the reference's taxonomic identity or biological annotation is correct.

**Important:** The guide must not contain COI primer sequences. Use the exact
verified coding interval and review the resulting alignment; a primer-free guide
is not a general guarantee that every primer residue is removed from all reads.
See [reference selection and preparation](coi_reference.md).

```bash
conda activate taxonomy
python code/postprocess_sintax.py --bold --kcalign \
  --kcalign-reference data/postprocessing/coi_refs/SATIP962-10_Nigrotipula_nigra.fasta \
  --taxon Tipulidae --folders barcode{01..21} \
  --analysis-dir data/sintax_analysis/mergedseq5_20261021 \
  --output-dir data/postprocessing/Tipulidae_bold_kcalign_run01 \
  --translation-table 5 --threads 6
```

The example requires a valid reference FASTA with a `>` header and a trusted,
correctly framed coding sequence. The module does not guess its frame or silently
trim bases. Divisibility by three and absence of internal stops alone do not prove
the biological reading frame; the user must verify its codon boundaries.

The alignment takes `combined_TAXON_nonumt_DATE.fasta`, the nucleotide consensus
sequences retained by the >=200 consecutive amino-acid screen. With
`--skip-translation`, it instead takes `combined_TAXON_selected_DATE.fasta`.
The corresponding prefix is `combined_TAXON_nonumt_DATE` or
`combined_TAXON_selected_DATE`. The mixed-mode command is equivalent to:

```bash
kc-align --mode mixed --reference trusted_reference.fasta \
  --sequences combined_TAXON_nonumt_DATE.fasta \
  --table 5 --keep --threads 6
```

**Installed-tool compatibility:** `--keep` is a switch, not a filename argument.
The tool writes `pre_align.fasta` and `codon_align.fasta` in its working directory.
Postprocessing isolates those files in a temporary directory, restores tracking headers,
and publishes `PREFIX.pre_align.aa.fasta` and `PREFIX.aln.fasta`. The installed
CLI omits the threads argument when calling mixed mode; `code/kcalign_cli.py`
forwards it correctly to the installed library without changing the conda package.
The local parallel adapter also fixes an installed rejection-handling bug that
otherwise crashes when a worker rejects a read. It runs the original serial
per-read checks in isolated worker directories, preserves rejected IDs, and uses
a persistent pool with at most `threads / 3` workers. Biological checks remain unchanged.
Mixed mode requires **1 thread or a positive multiple of 3**. The module passes
its `--threads` value unchanged and rejects incompatible values such as 4 or 8.
Kalign and MAFFT (the library's fallback aligner) must also be available in `taxonomy`.

The trusted reference guides alignment but is not added as an extra phylogenetic
tip. An independent input consensus with the same sequence remains eligible.
Mixed mode can reject reads with suspected frameshifts, early stops or excessive
ambiguity even after the earlier translation screen. `PREFIX.kcalign_selection.tsv`
records every input ID as `aligned` or `rejected_by_kcalign`; details are in
`PREFIX.kcalign.stdout.log` and `postprocessing.log`.
`PREFIX.kcalign_ids.tsv` maps temporary tool IDs to original headers.
At least four aligned consensus sequences are required for IQ-TREE. ClipKIT trims
the published **nucleotide** `PREFIX.aln.fasta` before IQ-TREE; the retained pre-alignment AA file
is an audit output. ClipKIT smart-gap column trimming follows alignment as described below.

MAFFT remains the default for general nucleotide alignment. kc-align is for coding
sequences; skipping the earlier SeqKit filter does not disable kc-align's own
coding-sequence checks. `--translation-table` controls its genetic code (default 5).

### ClipKIT smart-gap trimming

After either MAFFT or kc-align, Postprocessing automatically runs:

```bash
clipkit PREFIX.aln.fasta -m smart-gap -s nt -if fasta -of fasta \
  -o PREFIX.clipkit.aln.fasta
```

Smart-gap determines a gap threshold and removes gappy **alignment columns**.
It does not strip every gap independently from each sequence, and some gaps may
remain. The original alignment is preserved. Headers and the sequence set must
remain unchanged, and the trimmed alignment must be nonempty and rectangular.
IQ-TREE receives `PREFIX.clipkit.aln.fasta`. `clipkit.stdout.log` records the tool's
report; `clipkit.complete.json` records mode and input hash for resume validation.

This step is active for both aligners and with or without the translation screen.
No extra flag is required. On `--resume`, a completed untrimmed alignment can be
reused: missing or stale ClipKIT output is regenerated, and the dependent tree and
figure are rebuilt. Existing dependent outputs are archived. Failed partial trim
outputs are removed; extraction, clustering and alignment need not be repeated.

### Phylogeny and figure

For new runs, reference selection appends the source taxon label to the sequence
ID before length filtering, clustering and alignment: for example,
`BOLD_000000123__taxon__Tipula_maxima` or `AB123456.1__taxon__Tipula_maxima`.
Whitespace becomes underscores; punctuation unsafe for alignment/Newick tools
also becomes underscores. BOLD uses the species identification in its metadata; NCBI uses the
GenBank organism name. Original descriptions remain in the source files and audit
metadata, while selected reference headers contain only the annotated ID.
These names carry through both untrimmed and trimmed alignments into IQ-TREE.
Plotting accepts both original reference IDs and IDs with the `__taxon__` suffix,
using the original ID to look up source metadata. Local read headers are unchanged.
Resuming an older run reuses its existing sequence headers; this change does not
automatically rename completed runs.

IQ-TREE consumes the trimmed nucleotide `.clipkit.aln.fasta`, with `-m MFP`, seed 666
and the requested thread count. Bootstrap replicates are disabled; the pipeline
estimates a maximum-likelihood tree without bootstrap branch support. The
`--bootstrap-method` option is no longer accepted.
At least four retained consensus sequences are required. Outputs use `PREFIX`.

Baltic renders a midpoint-rooted divergence tree with branch lengths in
substitutions/site and its scale bar. Existing IQ-TREE support labels are displayed
on the corresponding splits, including slash-separated values. If `.treefile` has
no support and `.contree` does, the supported consensus topology is used. Supports
are preserved during rerooting; none are invented. Midpoint rooting does not
establish the biological root.

Tip labels are trusted taxon + read ID + folder + forward primer for local cluster seeds; GenBank
organism + accession.version for downloaded seeds; or available BOLD taxon + process
ID for BOLD seeds. GenBank metadata is fetched once and cached as
`reference_records.gb.xml`; concatenated Entrez XML batches are parsed locally.
BOLD metadata requires no network.

For the forward-primer field in local tip labels, `none` means tracking was
available but no forward primer was detected; `unknown` means the raw parent had
no matching tracking record. The label comes from the local cluster seed's
exported primer annotation. These values remain separate during clustering, as
described above; BOLD and GenBank tip labels do not include this field.

Figure outputs are `PREFIX.figure.pdf`, `.png`, `.svg`, `.midpoint.nwk`, `.tips.tsv` and
`.plot.json`. The tips table columns are `tree_tip_id`, `source`, `taxon_or_organism`,
`accession`, `readid`, `source_folder`, `display_label`; BOLD process IDs occupy
`accession` for BOLD rows. To plot independently use `code/plot_phylogeny.py --tree
PATH --analysis-dir PATH --reference-metadata PATH --output-prefix PATH`, with
GenBank XML or the BOLD metadata TSV as the metadata input.

### Resume and audit

Add `--resume` to the original command with the same output directory and settings.
Completed retained/selected sequences restart at alignment. A successfully completed
alignment restarts at ClipKIT, or at IQ-TREE if trimming is already complete,
without repeating extraction or reference retrieval.
Completion hashes distinguish finished alignments from interrupted outputs. MAFFT
partial files are removed on interruption; kc-align works in temporary files and
publishes its completion marker only after the nucleotide alignment validates.
Incomplete outputs are rebuilt on resume, with stale files archived under
`resume_archive/` and successful reference retrieval reused when needed.

Source, input database, taxon, folders, translation choice, aligner, genetic code
the trusted reference contents, and other biological settings must match the saved run. To change aligner or trusted reference while reusing validated retained sequences,
use `--resume --restart-alignment`. This archives the previous settings and rebuilds
alignment, tree and figure; subsequent resumes use the new settings. Other biological
settings must still match. Existing biological outputs are not rewritten by a code update.
To apply the new selection and clustering rules, rerun the analysis from extraction.

`parameters.json` stores settings, policy versions and executable paths;
`checkpoints.json` stores output hashes. `postprocessing.log` records commands,
stderr and counts. `postprocessing.failed` records failure; `postprocessing.complete`
is written only after phylogeny and plotting finish. DATE defaults to the run date
and may be specified with `--date YYYYMMDD`.

### Resuming runs with earlier bootstrap settings

The first `--resume` after disabling bootstrap reuses completed alignment and ClipKIT
outputs. The tree and figure are regenerated without bootstrap;
previous IQ-TREE outputs, including its checkpoint, are archived. A new
`phylogeny.complete.json` records the no-bootstrap setting and alignment hash.
Stop any still-running job for the same output folder before resuming; editing
this script does not stop an already-running IQ-TREE process. No branch-support
values are invented for the new unsupported tree.

BOLD references must have a named species identification. In full BOLD lineages,
the eighth comma-separated field is species (empty positions are preserved);
the optional ninth subspecies field is not used as the species label. Shortened
lineages are accepted only when their final nonempty entry is a named binomial.
Genus-only, missing and provisional identifications (`sp.`, `spp.`, `cf.`, `aff.`)
are excluded. Named binomials, optional subgenus annotations and trinomial names
are accepted. If no qualifying records match, extraction stops with an error.
