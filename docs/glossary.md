# Glossary

| Term | Meaning |
| --- | --- |
| **Alignment** | Sequences written one below the other, with gaps (`-`) inserted so that corresponding positions line up in columns. |
| **Barcode folder** | A folder such as `barcode01` written by MinKNOW: all reads carrying one ONT native barcode, i.e. one pool of PCR products. Here also called a *sample*. |
| **Bootstrap support** | A number on a tree branch (0–100) saying how consistently the data support that branch. Module 2's draft tree has none; add it for the final tree. |
| **BOLD** | Barcode of Life Data System, a database of DNA barcodes from identified specimens. |
| **Chimera** | An artificial read made of two molecules joined together. fastplong splits reads at internal adapters. |
| **ClipKIT** | Program that removes alignment columns that are mostly gaps. |
| **COI** | Cytochrome c oxidase subunit I, the mitochondrial gene used as the animal DNA barcode. *COI-5P* is the barcode region at its 5' end. |
| **Confidence** | SINTAX's support (0–1) for a name at one rank. Ranks with confidence ≥ 0.8 are *accepted*. |
| **Consensus sequence** | One sequence summarising a group of similar reads: at each position, the most common base. |
| **Contamination flag** | Warning that a taxon was found in ≥ 80% of the barcode folders. Not proof of contamination. |
| **Core (primer core)** | The part of a primer that binds the COI gene, without the tag. |
| **cutadapt** | Program that finds primer sequences in reads and cuts them off. |
| **Degenerate base** | A letter standing for several bases, used in primers: R = A/G, K = G/T, M = A/C, W = A/T, Y = C/T, N = any. |
| **Edit distance / error** | Number of substitutions, insertions and deletions needed to turn one sequence into another. |
| **fastplong** | Program that filters long reads by quality and length and removes adapters. |
| **FASTA / FASTQ** | Text formats for sequences. FASTA: a `>` header line then the sequence. FASTQ: also a quality score per base. `.gz` means compressed. |
| **Forward / reverse tag** | The 16-base label at the start of a tagged primer. The forward tag (`bc1001`–`bc1096`) identifies the specimen within a pool. |
| **IQ-TREE** | Program that estimates a maximum-likelihood phylogenetic tree. |
| **kc-align** | Aligner that aligns protein-coding DNA codon by codon, guided by one trusted sequence. |
| **MAFFT** | General-purpose sequence aligner. |
| **MIDORI** | Database of reference mitochondrial sequences with taxonomy, used by SINTAX. |
| **Midpoint rooting** | Placing the tree's root halfway along its longest path. A drawing choice, not the true ancestor. |
| **NUMT** | A nuclear copy of mitochondrial DNA. It can look like COI but evolves differently and misleads trees. |
| **ONT** | Oxford Nanopore Technologies, the sequencing platform. Its reads contain more insertion/deletion errors than short-read platforms. |
| **Reading frame / stop codon** | Protein-coding DNA is read in triplets (codons); a stop codon ends the protein. A frameshift or premature stop suggests an error, a pseudogene or another marker. |
| **SINTAX** | Classification method (in VSEARCH) that assigns each read a taxonomy with a confidence per rank. |
| **Tag block** | The first 21 bases of each tagged primer: the `GGTAG` spacer plus the 16-base tag. |
| **Taxid** | NCBI's number for a taxon, e.g. the `_12345` at the end of a MIDORI name. |
| **UGENE** | Desktop program to view and edit alignments. |
| **VSEARCH** | Program used here for classification (SINTAX) and for clustering reads into consensus sequences. |
