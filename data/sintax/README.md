# MIDORI reference database for Read processing

Run commands from the repository root with the `taxonomy` environment active.
The reference FASTA and generated VSEARCH `.udb` are large local files excluded
from Git; this README records how to obtain and place the reference.

## Download the matching release and format

Open the official [MIDORI download page](https://www.reference-midori.info/download.html).
It links to the current release on Zenodo and a separate archive of previous
releases. Use **Previous releases**, then select **GB272 (June 2026)** to reproduce
this repository's default configuration. For another release, keep its real
release number in the filename and change the configuration accordingly.

Select **SINTAX**, **UNIQ**, **NUC_SP**, **CO1**. NUC_SP contains nucleotide
records including ambiguous species names; CO1 is the COI marker used here.
Download the archive containing that combination and extract the matching FASTA:

```text
MIDORI2_UNIQ_NUC_SP_GB272_CO1_SINTAX.fasta
```

Keep the taxonomy-bearing FASTA headers intact. BLAST, RAW, amino-acid and other
classifier formats are not substitutes for the SINTAX file. Record the release,
download URL and any supplied checksum when preparing a new dataset. Follow the
provider's citation instructions.

## Folder and decompression

Create the destination:

```bash
mkdir -p data/sintax
```

If the download is an archive containing multiple files, unpack it first and
copy only the selected CO1 SINTAX FASTA here. The pipeline does not unpack ZIP or
TAR archives. Once any outer archive is extracted, either layout is accepted:

```text
data/sintax/MIDORI2_UNIQ_NUC_SP_GB272_CO1_SINTAX.fasta
```

or, if the selected FASTA itself is gzipped:

```text
data/sintax/MIDORI2_UNIQ_NUC_SP_GB272_CO1_SINTAX.fasta.gz
```

Manual gzip decompression is optional. When classification is needed and the
`.udb` is absent, `code/run_metabarcoding_pipeline.sh` decompresses the matching
`.fasta.gz` if the plain FASTA is absent, then builds the UDB automatically. It
keeps the gzip original. Allow disk space for both the plain FASTA and the UDB.
For manual decompression while keeping the gzip file:

```bash
gzip -dk data/sintax/MIDORI2_UNIQ_NUC_SP_GB272_CO1_SINTAX.fasta.gz
```

The settings in `config/pipeline.env` refer to the **uncompressed** FASTA path,
even when initially supplying only its `.gz` sibling:

```bash
SINTAX_REFERENCE=data/sintax/MIDORI2_UNIQ_NUC_SP_GB272_CO1_SINTAX.fasta
SINTAX_DB=data/sintax/MIDORI2_UNIQ_NUC_SP_GB272_CO1_SINTAX.udb
```

An existing nonempty UDB is reused without checking whether its FASTA changed.
For a new database release, use a new matching FASTA/UDB pair and distinct
classification outputs (`SINTAX_OUTPUT_ROOT` and `SINTAX_SUFFIX`); update the
post-classification settings to read those outputs. Do not rename a new release
as GB272 or reuse an old UDB with it.

SINTAX sample tables under `data/sintax/DATASET/barcode*/` and the UDB are generated
outputs, not additional downloads. See the [Read processing guide](../../docs/pipeline_from_fastplong_onwards.md).

Classification outputs now default to a separate subdirectory per dataset. For
older outputs directly under `data/sintax/barcode*/`, explicitly select
`--sintax-root data/sintax` in the read-processing/report launcher.
