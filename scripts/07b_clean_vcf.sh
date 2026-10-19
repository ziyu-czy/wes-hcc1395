#!/bin/bash
# Step 7b: remove the normal-LOH regions (6p, 16q, chrX) from the PASS VCF. Conda env: wes.
# To inspect the leaked sites first:
#   bcftools query -R vep/germ_pos.tsv -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE:%AD:%AF]\n' vcf/HCC1395.pass.vcf.gz | head -20 | column -t
set -euo pipefail
bcftools view -T ^vcf/normal_loh.bed vcf/HCC1395.pass.vcf.gz -Oz -o vcf/HCC1395.clean.vcf.gz
bcftools index -t vcf/HCC1395.clean.vcf.gz
bcftools view -H vcf/HCC1395.clean.vcf.gz | wc -l   # 1,270
