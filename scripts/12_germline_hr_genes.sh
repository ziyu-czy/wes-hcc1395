#!/bin/bash
# Step 12: germline variants in homologous-recombination genes (BRCA1, BRCA2, PALB2, RAD51C).
# Mutect2 removes germline variants by design, so an inherited first hit would be invisible in steps 4-9.
# Here both samples are genotyped together with bcftools. Conda env: wes.
set -euo pipefail
REF=ref/Homo_sapiens_assembly38.fasta
REGIONS=chr17:43044000-43126000,chr13:32315000-32400300,chr16:23603000-23642000,chr17:58692000-58736000
#        BRCA1                   BRCA2                   PALB2                   RAD51C

bcftools mpileup -f $REF -a AD,DP -q 20 -Q 20 -r $REGIONS \
  bam/HCC1395BL.dedup.bam bam/HCC1395.dedup.bam \
  | bcftools call -mv -Oz -o vcf/hr_genes_germline.vcf.gz

# keep sites with depth >= 10 in the normal (sample 0)
bcftools view -i 'FMT/DP[0]>=10' vcf/hr_genes_germline.vcf.gz -Oz -o vcf/hr_genes_germline.dp10.vcf.gz
bcftools view -H vcf/hr_genes_germline.dp10.vcf.gz | wc -l   # 50

# genotype table for joining with the VEP output (normal first, then tumor)
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%AD]\n' \
  vcf/hr_genes_germline.dp10.vcf.gz > vep/hr_genes_gt.tsv

# Next: annotate vcf/hr_genes_germline.dp10.vcf.gz with Ensembl VEP (web), selecting HGVS, AlphaMissense,
# ClinPred, SpliceAI, gnomAD Exomes and "Phenotype & disease associations" (provides the ClinVar_* columns).
# Save as vep/hr_genes_germline_VEP.tsv.gz, then run scripts/12b_germline_hr_genes.R
