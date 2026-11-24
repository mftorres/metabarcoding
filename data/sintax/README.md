# MIDORI reference database for Module 1

Module 1 classifies reads with VSEARCH SINTAX against the MIDORI COI database.
The database is large and stays on your computer; it is not stored in Git.

## Download the right file

1. Open the [MIDORI download page](https://www.reference-midori.info/download.html).
2. To reproduce this project's settings, choose **Previous releases**, then
   **GB272 (June 2026)**. (For a newer release, keep its real release number in the
   file name and update `config/pipeline.env`, see below.)
3. Select the combination **SINTAX · UNIQ · NUC_SP · CO1** and download it.
   - *SINTAX* is the file format VSEARCH needs; other formats (BLAST, RAW…) do not work.
   - *NUC_SP* contains nucleotide sequences, including records with ambiguous species names.
   - *CO1* is the COI marker.
4. If the download is a ZIP or TAR archive, unpack it and keep only
   `MIDORI2_UNIQ_NUC_SP_GB272_CO1_SINTAX.fasta` (or its `.fasta.gz`).

Write down the release, download date and link. Cite MIDORI as its website asks.

## Put it in place

```bash
cd ~/metabarcoding
mkdir -p data/sintax
mv "/mnt/c/Users/USERNAME/Downloads/MIDORI2_UNIQ_NUC_SP_GB272_CO1_SINTAX.fasta.gz" data/sintax/
```

Either `data/sintax/MIDORI2_UNIQ_NUC_SP_GB272_CO1_SINTAX.fasta` or the same name with
`.gz` works. The first time classification runs, the pipeline unpacks a `.gz` copy (it
keeps the original) and builds a VSEARCH index file, `…_SINTAX.udb`. This takes a while
and needs free disk space for both: the GB272 index is about 10 GB.

## Using another release

Change these three settings in `config/pipeline.env` together, so new results are not
mixed with old ones:

```bash
SINTAX_REFERENCE=data/sintax/MIDORI2_UNIQ_NUC_SP_GB999_CO1_SINTAX.fasta
SINTAX_DB=data/sintax/MIDORI2_UNIQ_NUC_SP_GB999_CO1_SINTAX.udb
SINTAX_SUFFIX=sintaxGB999
```

An existing `.udb` is reused without checking which FASTA it came from, so never
rename a new release to an old name.

## Optional: an Arthropoda-only database

Keeping only arthropod records makes the database smaller and SINTAX faster. Done for
this project on 23 Nov 2026, in the `taxonomy` environment, from `data/sintax/`:

```bash
# Keep records whose phylum field is exactly Arthropoda (NCBI taxid 6656).
seqkit grep -n -r -p ',p:Arthropoda_6656,' \
  MIDORI2_UNIQ_NUC_SP_GB272_CO1_SINTAX.fasta > MIDORI2_UNIQ_NUC_SP_GB272_CO1_ARTH_SINTAX.fasta
grep -c '^>' MIDORI2_UNIQ_NUC_SP_GB272_CO1_ARTH_SINTAX.fasta   # GB272: 2,432,079 of 3,191,213

# Build the VSEARCH index for the subset.
vsearch --makeudb_usearch MIDORI2_UNIQ_NUC_SP_GB272_CO1_ARTH_SINTAX.fasta \
  --output MIDORI2_UNIQ_NUC_SP_GB272_CO1_ARTH_SINTAX.udb \
  --notrunclabels
```

- **Match the phylum field, not the word.** MIDORI writes the phylum as
  `p:Arthropoda_6656`. A pattern such as `\bArthropoda\b` does not match it, because
  regular expressions count `_` as part of a word. That pattern keeps only the ~1,000
  unidentified `Arthropoda sp.` records, and the index looks fine but is useless.
  Always check the record count with `grep -c '^>'`, as above.
- **Expect about a quarter less work.** Arthropods are 76% of the GB272 COI records,
  so SINTAX time drops by roughly that share, not several-fold.
- **Non-arthropod reads get forced into arthropods.** SINTAX can only choose among the
  names in the database. Reads from fungi, plants, vertebrates or bacteria are then
  given an arthropod name, sometimes with high confidence at high ranks. Use the
  subset only when such reads are of no interest. Do not use it when you need to
  report or remove them.
- **Keep the results apart.** Set the three settings together in `config/pipeline.env`,
  for example `SINTAX_SUFFIX=sintaxARTH`, and pass `--suffix sintaxARTH` to
  `sintax_analysis.py`. Results made with the full database are named `…_sintaxSP.tsv`.

Classification results are written to `data/sintax/DATASET/barcodeNN/barcodeNN_sintaxSP.tsv`
(see [Module 1](../../docs/read_processing.md)).
