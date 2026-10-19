# Tumor–normal whole-exome analysis of HCC1395: from somatic calls to a BRCA1-deficient genome

An end-to-end analysis of a public tumor–normal exome pair: the triple-negative breast cancer cell line **HCC1395** and the matched lymphoblastoid line **HCC1395BL** from the same donor, both sequenced by the [SEQC2 consortium](https://doi.org/10.1038/s41587-021-00993-6). The pipeline goes from raw FASTQ to somatic mutations, benchmarking, mutational signatures, allele-specific copy number and germline variants. Along the way I followed up on every anomaly I hit, and those follow-ups produced most of the biology in this report.

## Key findings

- **Somatic SNV calls are accurate inside the SEQC2 high-confidence regions** (precision 0.969, recall 0.898). Indel calls are much weaker (precision 0.473).
- **The normal sample is not a clean reference.** HCC1395BL has lost one haplotype on 6p, 16q and chrX. As a result, about 340 germline SNPs were called as somatic. The SEQC2 benchmark cannot detect this because those regions sit outside its high-confidence set.
- **Driver mutations:** *TP53* p.Arg175His, with loss of the wild-type allele through LOH of chromosome 17. There is also a somatic *BRCA2* p.Glu1593Ter, but it is monoallelic: the wild-type allele is still present.
- **Mutational processes:** APOBEC (SBS2 + SBS13) accounts for about one third of SNVs. SBS3, the HR-deficiency signature, cannot be distinguished from the flat clock-like signatures with about 1,100 SNVs.
- **Genome structure:** the tumor is near-triploid (mean ploidy about 2.8) with extensive LOH, a pattern consistent with whole-genome doubling. The HRD-LOH score is elevated (≥ 10).
- **Cause of HR deficiency:** a **germline pathogenic *BRCA1* c.5251C>T p.Arg1751Ter** (ClinVar: Pathogenic, expert panel). The same chr17 LOH removed the wild-type *BRCA1* and the wild-type *TP53*.
- **Copy-number drivers:** homozygous deletion of ***CDKN2A***; partial homozygous deletion of ***PTEN***; gains of *MYC* and *CCND1*. *ERBB2* is not amplified, consistent with a HER2-negative tumor.

## Data

| | Tumor | Normal |
|---|---|---|
| Sample | HCC1395 (breast carcinoma, TNBC) | HCC1395BL (EBV-transformed B-lymphoblastoid, same donor) |
| SEQC2 library / SRA run | WES_IL_T_1 / SRR7890883 | WES_IL_N_1 / SRR7890874 |

- Illumina HiSeq 4000, 2 × 125 bp. FASTQs were downloaded from ENA and md5-verified.
- Capture kit: Agilent SureSelect Human All Exon V6 + UTR (S07604624). The hg38 target BED comes from the SEQC2 FTP.
- Reference: GATK hg38 bundle (`Homo_sapiens_assembly38.fasta` + ALT contigs). Mutect2 resources: `af-only-gnomad.hg38`, `1000g_pon.hg38`, `small_exac_common_3.hg38`.
- Truth set: SEQC2 high-confidence somatic SNV/indel v1.2.1 and High-Confidence Regions v1.2.

## Pipeline

| Step | What it answers | Tools | Script |
|---|---|---|---|
| 1 | Read QC and alignment | FastQC, bwa-mem (ALT-aware), samtools | `01_align.sh` |
| 2 | Duplicate marking | GATK MarkDuplicates | `02_markdup.sh` |
| 3 | Exome target intervals | bedtools | `03_targets.sh` |
| 4–5 | Somatic SNV/indel calling and filtering | GATK Mutect2, LearnReadOrientationModel, CalculateContamination, FilterMutectCalls | `04_mutect2.sh`, `05_filter.sh` |
| 6 | How accurate are the calls? | bcftools norm/isec vs SEQC2 truth | `06_benchmark.sh` |
| 7 | Annotation and a germline-leakage check | Ensembl VEP (web), R | `07_annotate_and_normal_loh.R`, `07b_clean_vcf.sh` |
| 8 | Which drivers are mutated? | R | `08_drivers.R` |
| 9 | Which mutational processes acted? | MutationalPatterns (context from local FASTA via Rsamtools) | `09_signatures.R` |
| 10 | Copy-number ratios and segments | CNVkit | `10_cnvkit.sh` |
| 11 | Ploidy, BAF/LOH, HRD-LOH | Python (pandas, matplotlib) | `11_cnv_baf_hrd.py` |
| 12 | Is there an inherited first hit in HR genes? | bcftools mpileup/call, VEP (ClinVar) | `12_germline_hr_genes.sh`, `12b_germline_hr_genes.R` |
| 13 | Amplified or deleted cancer genes | Python, samtools bedcov | `13_cn_drivers.py`, `13b_check_deletions.sh` |

BQSR was skipped on purpose. Base qualities on this HiSeq 4000 data are already well calibrated, its effect on somatic SNV calling is small, and the time it saves on a laptop matters.

## Results

### 1. Alignment and somatic calling

Both samples map at more than 99.7%. The tumor has fewer properly paired reads than the normal (93.5% vs 98.9%) and more read pairs spanning two chromosomes (4.8% vs 0.14%). The shared hotspots come from reference repeats, while the tumor's excess is spread thinly across the genome. This points to chimeric fragments from library preparation rather than real translocations. Both libraries have about 43% duplicates, which reflects sequencing saturation. Mutect2 called 10,219 candidates, and **1,750** passed all filters. Estimated contamination was 0.12%.

### 2. Benchmark against SEQC2

Evaluated within the exome targets ∩ SEQC2 high-confidence regions (80.7 Mb):

| | TP | FP | FN | Precision | Recall |
|---|---|---|---|---|---|
| SNV | 1,041 | 33 | 118 | 0.969 | 0.898 |
| Indel | 43 | 48 | 7 | 0.473 | 0.860 |

Most indel false positives are long, weakly supported events (10–90 bp). They are likely assembly artifacts from soft-clipped or chimeric reads. I did **not** tune filters against this truth set: using the same truth set to tune and to evaluate would overstate performance.

### 3. A hidden problem: LOH in the normal cell line

True somatic mutations should almost never be common population variants. Yet **276 PASS calls had gnomAD allele frequency > 0.001**, some as high as 0.75, and **272 of them clustered in three blocks**: 6p (1.6–53 Mb), all of 16q, and all of chrX. At these sites the normal showed **zero alt reads at 27–205× depth**, while the tumor showed a consistent VAF of about 0.65.

The explanation is that HCC1395BL, an EBV-transformed line grown in culture, **lost one haplotype in these regions**. The matched normal therefore confidently "proved" that the donor did not carry these SNPs, and Mutect2 called them as somatic. gnomAD only adjusts the prior, so it cannot override direct evidence from a deeply covered normal.

- Inside the three regions, 67% of calls were in gnomAD, compared with 20% elsewhere. I estimate that about 340 of the 480 calls there are germline. I excluded these regions, which leaves a **clean set of 1,270 variants**.
- Only 1 of the 276 leaked sites is in the SEQC2 truth set, and only 2 fall inside the evaluation regions. **SEQC2 excluded these regions from its high-confidence set, so the benchmark is blind to this problem.** A good benchmark score does not guarantee a clean genome-wide call set.
- In step 11 the BAF data later confirmed the same thing independently: chrX has 581 informative SNP sites, but not one of them is heterozygous in the normal.

### 4. Driver mutations

Of the 364 protein-altering variants in the clean set, two hit established breast-cancer drivers:

| Gene | Change | Tumor AD (ref, alt) | Normal AD | Interpretation |
|---|---|---|---|---|
| *TP53* | p.Arg175His (AlphaMissense 0.986) | 0, 213 | 323, 0 | hotspot; VAF ≈ 1 → wild-type allele lost |
| *BRCA2* | p.Glu1593Ter | 64, 46 (VAF 0.42) | 86, 0 | truncating but **monoallelic**; no second small variant, even among filtered candidates |

### 5. Mutational signatures

I built the 96-class SBS profile from the 1,140 clean SNVs. The trinucleotide context was read directly from the local reference FASTA, which made the 3 GB BSgenome package unnecessary.

- C>G (31%) and C>T (35%) dominate, with sharp peaks at **TCA/TCT**. This is the APOBEC signature (SBS2 + SBS13), and it is stable at about one third of SNVs in every fit. There is no C>A-dominated oxidative artifact, which is a sign that the calls are good.
- A strict refit against all 60 COSMIC signatures over-fits. It assigns UV (SBS7a), platinum (SBS31) and haloalkane (SBS42) signatures to a breast cancer cell line.
- **Testing SBS3 directly:** I fitted the breast-relevant signatures with and without SBS3. Without it, SBS40/SBS5 absorb its whole share and the fit barely changes (cosine 0.962 → 0.955). Flat signatures cannot be told apart with this few mutations, so **signatures can neither confirm nor exclude HR deficiency here**. Answering the question needed a different line of evidence.

### 6. Copy number and allelic imbalance

- **Ploidy.** CNVkit segment levels pile up at log2 ≈ −0.5, +0.08 and +0.5. The spacing matches 2, 3 and 4 copies, so log2 = 0 corresponds to the tumor's average ploidy of about 2.8, i.e. near-triploid. The full spacing also shows purity close to 100%. CNVkit's own integer calls, which assume a diploid genome, were therefore not used.
- **Noise.** The depth ratio is noisy (MAD 0.33, bivar 0.35), and target and off-target bins are about equally noisy. The likely cause is the quality difference between the two libraries combined with having only one normal and no panel of normals.
- **BAF.** BAF was measured at 4,753 SNPs that are heterozygous in the normal. Because it is a ratio within one sample, it is immune to that library noise. It shows **extensive LOH**: whole chromosomes 3, 5, 8, 10 and 17, plus large parts of others. Many LOH regions are 2:0, meaning two copies of the same haplotype. A loss followed by whole-genome doubling explains this pattern.
- **chr17.** The whole chromosome is LOH, and 17p carries about 3 copies, all of them with *TP53* R175H. The mutation therefore arose before the region was copied.
- **BRCA2 region.** Both haplotypes are present (≈ 4:1), so there is no LOH at *BRCA2*.
- **HRD-LOH** (LOH regions > 15 Mb that are not whole chromosomes) = **20**. A conservative re-count that merges fragmented regions and drops borderline ones gives 10–14. Both are at or above the usual threshold of 10. One caveat: near-triploid genomes tend to score higher.

### 7. Why is HR deficient? A germline *BRCA1* first hit

The HR scars needed a cause, and the monoallelic *BRCA2* could not explain them: the scars are present in every tumor cell, while wild-type *BRCA2* is still present. Mutect2 discards germline variants by design, so I genotyped both samples at *BRCA1*, *BRCA2*, *PALB2* and *RAD51C* with bcftools and annotated the results with ClinVar.

| Gene | Variant | ClinVar | gnomAD AF | Normal AD | Tumor AD |
|---|---|---|---|---|---|
| ***BRCA1*** | **c.5251C>T p.Arg1751Ter** (stop gained) | **Pathogenic** (reviewed by expert panel) | 2.7 × 10⁻⁶ | 86, 71 (het) | **0, 67** |

Every other variant that was heterozygous in the normal is ClinVar Benign. The picture is the classic two-hit model: **an inherited pathogenic *BRCA1* allele plus somatic loss of the wild-type allele through chr17 LOH**. That single LOH event also completed the inactivation of *TP53*.

### 8. Copy-number drivers

| Gene | Estimated copies | Interpretation |
|---|---|---|
| *CDKN2A* | 0 | **homozygous deletion**: every exon ≈ 0× in the tumor vs 47–284× in the normal |
| *PTEN* | 0 for the 3′ part | **partial homozygous deletion**: 5′ exons retained, ≈ 0× from about chr10:87.957 Mb to the gene end, on an LOH background |
| *MYC* | ≈ 6.4 | gain (about 2× mean ploidy) |
| *CCND1* | ≈ 5.5 | gain |
| *ERBB2* | ≈ 2 (LOH) | not amplified → consistent with HER2-negative TNBC |
| *FGFR1*, *EGFR*, *PIK3CA*, *RB1* | 2–4 | near mean ploidy |

The *PTEN* deletion is shorter than CNVkit's segmentation resolution. Its segment looked like 2 copies, and only the gene's own bins and the exon-level read depth in the BAMs revealed it.

### Final picture

| Category | Alteration | How |
|---|---|---|
| HR deficiency | *BRCA1* p.Arg1751Ter (germline) | + chr17 LOH → biallelic |
| Tumor suppressor | *TP53* p.Arg175His | + the same chr17 LOH |
| Tumor suppressor | *CDKN2A* | homozygous deletion |
| Tumor suppressor | *PTEN* | LOH + partial homozygous deletion |
| Oncogene | *MYC*, *CCND1* | copy gain |
| Passenger / later event | *BRCA2* p.Glu1593Ter | somatic, monoallelic |
| Mutational process | APOBEC (SBS2/SBS13) | ≈ 1/3 of SNVs |
| Genome | near-triploid, extensive LOH, HRD-LOH ≥ 10 | consistent with HR deficiency and whole-genome doubling |

Taken together, this is the genome of a typical ***BRCA1*-associated triple-negative breast cancer**. Independent vendor data also list HCC1395 as homozygous for *TP53* p.Arg175His.

## Limitations

- **Single normal, no panel of normals.** The depth-based copy number is noisy, and segment boundaries and small events are uncertain. The ploidy estimate rests on segment spacing plus BAF rather than on a dedicated purity/ploidy model such as FACETS or ASCAT.
- **The normal is a cell line with its own LOH.** Calls in 6p, 16q and chrX were removed, so somatic events in those regions are missed. BAF is also uninformative there.
- **Indel calling** is weak (precision 0.47), and I did not optimise it.
- **Signatures** come from about 1,100 exome SNVs, and the COSMIC signatures were not renormalised for exome trinucleotide composition.
- **HRD-LOH** is a simplified re-implementation on noisy segments. I did not compute TAI or LST, so this is not a validated HRD score.
- **Methylation and structural variants** are invisible to WES SNV/CNV analysis. For example, *BRCA1* promoter methylation or the breakpoints of the *PTEN* deletion cannot be seen here.
- **Annotation** used the Ensembl VEP web interface, so it is not fully scripted.

## Repository structure

```
README.md
docs/file_formats.md     # what each file type is and where it appears in this pipeline
scripts/                 # numbered steps, run from the project root
results/tables/          # small result tables (benchmark, drivers, LOH regions, exon depths)
results/figures/         # 96-class profile, CNVkit scatter, segment histogram, BAF plots
```

Raw data (FASTQ, BAM, reference, VCF) are not included; see *Data* for their sources.

## Reproducing

1. Create the conda environments: `wes` (fastqc, bwa, samtools, bcftools, bedtools, gatk4) and one with `cnvkit`. Install the R packages `readr`, `dplyr`, `MutationalPatterns`, `VariantAnnotation` and `Rsamtools`, and the Python packages `pandas` and `matplotlib`.
2. Download the FASTQs (ENA: SRR7890883, SRR7890874), the GATK hg38 bundle and Mutect2 resources, the SEQC2 truth set and target BED, and the GENCODE 50 GTF (`ref/gencode.v50.annotation.gtf.gz`).
3. Run `scripts/01_*` → `13b_*` in order from the project root. Steps 7 and 12 include a manual Ensembl VEP web run; the settings are described in those scripts.

## Software

bwa 0.7.19 · samtools 1.24 · bcftools · bedtools 2.31.1 · GATK 4.6.2 · CNVkit 0.9.14 · Ensembl VEP (web, GENCODE 50) · R: MutationalPatterns, VariantAnnotation, Rsamtools, dplyr · Python: pandas, matplotlib

## Acknowledgements

Data: Sequencing Quality Control Phase 2 (SEQC2) Somatic Mutation Working Group. Fang LT *et al.*, *Nature Biotechnology* 39, 1151–1160 (2021). https://doi.org/10.1038/s41587-021-00993-6
