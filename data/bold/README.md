# BOLD reference database for Postprocessing

Run commands from the repository root with the `taxonomy` environment active.
The database stays local and is excluded from Git; commit this README instead.
Postprocessing needs it only for `--bold`; `--download` obtains references from NCBI.

## Download a FASTA snapshot

Open the official [BOLD Data Packages page](https://boldsystems.org/data/data-packages/),
choose a reference-library snapshot, and download its **FASTA sequence** file.
Keep the snapshot date, source link and accompanying citation/checksum information.
This checkout uses the local filename `BOLD_Public.16-Oct-2026.fasta.gz`; that is
an example snapshot name, not a permanent URL or a requirement to rename newer
releases. TSV/JSON metadata downloads cannot replace the sequence FASTA.

Postprocessing expects headers with this structure; preserve the full headers:

```text
>processID|marker|country|comma-separated lineage
```

## Placement and compression

```bash
mkdir -p data/bold
```

Place the prepared sequence database here, for example:

```text
data/bold/BOLD_Public.16-Oct-2026.fasta.gz
```

**Keep the FASTA gzipped:** Postprocessing streams gzip directly, so decompression is
not needed. Directory discovery accepts `.fasta.gz`, `.fa.gz` and `.fna.gz`.
The reader also expects gzip when an explicit `--bold-db` file is supplied.
If the download is plain FASTA, compress it with `gzip -k downloaded.fasta`, or
use the SeqKit preparation command below. If an outer ZIP/TAR package was
supplied, extract the sequence FASTA from that package first.

## Remove records without sequences

Some downloaded records may have a header but no sequence. A matching empty
record stops BOLD extraction. [SeqKit](https://bioinf.shenwei.me/seqkit/usage/#seq)
can discard zero-length records with `seq --min-len 1` and write gzip directly.
Keep the original download **outside `data/bold/`**, for example under
`/path/to/downloads/`, then prepare the file that Postprocessing will read:

```bash
seqkit seq --min-len 1 \
  /path/to/downloads/BOLD_Public.16-Oct-2026.fasta.gz \
  -o data/bold/BOLD_Public.16-Oct-2026.fasta.gz
seqkit stats data/bold/BOLD_Public.16-Oct-2026.fasta.gz
```

Replace the example date and input path with your downloaded snapshot. Plain
FASTA input works too. The input and output must be different files. This removes
empty records, preserves nonempty sequences and full headers, and does not select
a marker or impose Postprocessing's later 500–2,000 nt length filter.

## Select the database

By default, `--bold` scans every supported gzip FASTA in `data/bold/`. Keep only
the intended prepared snapshot(s) there; keeping raw and cleaned copies together
would select records twice. To choose one file explicitly:

```bash
python code/postprocess_sintax.py --bold \
  --bold-db data/bold/BOLD_Public.16-Oct-2026.fasta.gz \
  --taxon Tipulidae --folders barcode{01..21} \
  --analysis-dir data/sintax_analysis/mergedseq5_20261021 \
  --output-dir data/postprocessing/Tipulidae_bold_new_run --threads 8
```

Use a taxon name, not a numeric taxid. Header matching includes all matching
markers, so use a marker-specific snapshot if needed. The pipeline generates
selected references and metadata inside each run directory; those outputs do not
belong in this database folder. See the [Postprocessing guide](../../docs/post_processing.md).
