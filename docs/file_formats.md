# File formats in this pipeline, in the order you meet them

The whole pipeline is a chain of translations:

> **reads** (FASTQ) → **where each read sits on the genome** (BAM) → **where the genome differs from the reference** (VCF) → **what those differences mean** (annotation tables) → **how much DNA each region has** (CNVkit files)

Three "map" files are used along the way: the reference genome itself (FASTA), the list of regions we care about (BED), and the list of where genes and exons are (GTF).

---

## Reference files: the map

### FASTA (`.fasta`, `.fa`) – the reference genome sequence
- **What:** the reference DNA sequence, one entry per chromosome.
  ```
  >chr17
  GATCAATGAGGTGGACACCAGAGGCGGGGACTTGTAAATAACACTGGGCTGT...
  ```
- **Here:** `ref/Homo_sapiens_assembly38.fasta` (GATK hg38). Every coordinate in every other file refers to this sequence, so all files must use the same genome build.
- **Companion files:** `.fai` (index, so tools can jump to any position), `.dict` (chromosome names and lengths, needed by GATK), and the bwa index (`.bwt`, `.sa`, ...) used for alignment.
- **Tools:** `samtools faidx`; `Rsamtools::scanFa` (step 9 reads trinucleotide context from it).

### BED (`.bed`) – a list of genomic intervals
- **What:** one interval per line: chromosome, start, end (plus optional columns such as a name).
  ```
  chr1    69036    71585
  ```
- **Coordinates:** **0-based, half-open**. The first base of a chromosome is 0 and the end is not included. This differs from VCF and GTF, which are 1-based.
- **Here:** exome capture targets (`ref/exome_targets.hg38.bed`), SEQC2 high-confidence regions, the excluded normal-LOH regions (`vcf/normal_loh.bed`), and exon intervals for depth checks.
- **Tools:** `bedtools` (intersect, merge, sort); passed to other tools to restrict analysis (`-L`, `-R`, `-T`).

### GTF (`.gtf.gz`) – gene annotation
- **What:** the coordinates of every gene, transcript and exon, with attributes in the last column.
  ```
  chr10  HAVANA  exon  87863114  87864548  .  +  .  gene_name "PTEN"; exon_number 1; tag "MANE_Select"; ...
  ```
- **Coordinates:** 1-based, inclusive. Subtract 1 from the start when converting to BED.
- **Here:** `ref/gencode.v50.annotation.gtf.gz`, the same GENCODE release that VEP used. It provides gene spans and exon numbers (steps 13 and 13b).
- **Tools:** `grep`/`awk` for quick extraction; `rtracklayer`/`GenomicFeatures` in R.

---

## Raw data: what the sequencer produced

### FASTQ (`.fastq.gz`) – sequencing reads
- **What:** 4 lines per read: name, sequence, `+`, and per-base quality. Qualities are ASCII characters (Phred + 33).
  ```
  @SRR7890883.1
  NTGCTCAGGTTCAGCCTGGAGAAGGTGACCT...
  +
  #AAFFJJJJJJJJJJJJJJJJJJJJJJJJJJ...
  ```
- **Here:** paired-end, so `_1` and `_2` files per sample, 2 × 125 bp.
- **Tools:** FastQC (quality), bwa (alignment).

---

## Alignment: where each read sits

### SAM / BAM (`.bam`) + index (`.bai`)
- **What:** every read with its alignment: chromosome, position, mapping quality (MAPQ), CIGAR (how the read matches: `100M2I23M`), FLAG (paired? reverse strand? duplicate?), mate position, and tags (read group, NM, MD, ...). BAM is the compressed binary form of SAM; the `.bai` index allows random access.
- **Here:** `bam/*.sorted.bam` after bwa, `bam/*.dedup.bam` after duplicate marking. This is the evidence everything else is computed from.
- **Tools:** `samtools` (view, flagstat, stats, index, bedcov), GATK, CNVkit.
- **Key idea:** a BAM can show what a VCF cannot. A homozygous deletion leaves *no reads*, so step 13b counts depth in the BAM instead of looking for variants.

---

## Variants: where the genome differs

### VCF (`.vcf.gz`) + index (`.tbi`)
- **What:** one line per variant site: chromosome, position (1-based), REF, ALT, FILTER, INFO, then per-sample FORMAT fields such as genotype (GT), allele depths (AD, e.g. `64,46` = 64 ref and 46 alt reads) and allele fraction (AF).
  ```
  chr13  32339132  .  G  T  .  PASS  ...  GT:AD:AF  0/1:64,46:0.42  0/0:86,0:0.011
  ```
- **Indels:** POS points to the base *before* the change. VEP reports the first changed base, which is one position later (handled in step 12b).
- **Here:** Mutect2 output (unfiltered → filtered → PASS → clean). These contain somatic candidates only, because germline variants are removed by design. The bcftools germline calls in HR genes (step 12) contain both.
- **Tools:** `bcftools` (view, query, norm, isec, index, mpileup, call); `VariantAnnotation` in R.

### GATK pileup table (`.pileups.table`)
- **What:** ref and alt read counts in one sample at a fixed list of common SNPs.
- **Here:** made for contamination estimation (step 5) and re-used to compute **BAF** (step 11): find SNPs heterozygous in the normal, then read their allele fraction in the tumor.

---

## Annotation: what the differences mean

### VEP output (`.tsv.gz`)
- **What:** a table with **one row per variant × transcript**, which is why 1,750 variants gave 41,040 rows. Columns give the gene (`SYMBOL`), consequence and `IMPACT`, protein change (`HGVSp`, e.g. `p.Arg175His`), MANE/canonical flags, predictors (AlphaMissense, ClinPred, SpliceAI), gnomAD frequency and ClinVar significance.
- **Tips:** pick one row per variant (MANE Select first). Read all columns as text, then convert numbers explicitly. Mutect2 VCFs have no IDs, so build a key from `Location` + `Allele`.
- **Tools:** R (`readr`, `dplyr`). Use `dplyr::` prefixes when Bioconductor packages are loaded.

---

## Copy number: how much DNA each region has (CNVkit)

| File | Content |
|---|---|
| `*.target.bed`, `*.antitarget.bed` | bins on the capture targets, and larger bins in the off-target space between them |
| `*coverage.cnn` | read depth per bin, per sample |
| `reference.cnn` | expected depth per bin from the normal, with GC and other bias corrections |
| `.cnr` | **log2(tumor / reference) for every bin**: raw and noisy |
| `.cns` | **segments**: neighbouring bins with similar log2 merged into one value |
| `.call.cns` | integer copy number assuming a diploid, pure sample (not used here; the tumor is near-triploid) |

- **Key idea:** log2 = 0 means "equal to this tumor's average", not "2 copies". BAF (from the pileup tables) is needed to anchor absolute copy number and to see LOH.
- **Tools:** `cnvkit.py` on the command line; pandas (or `cnvlib`) to analyse the tables.

---

## Index files at a glance

| Data file | Index | Made by |
|---|---|---|
| FASTA | `.fai`, `.dict`, bwa index | `samtools faidx`, `gatk CreateSequenceDictionary`, `bwa index` |
| BAM | `.bai` | `samtools index` |
| VCF.gz | `.tbi` | `bcftools index -t` |

An index lets a tool jump straight to a region (`-r chr17:7675088`) instead of reading the whole file.
