#!/bin/bash
# Step 5: model orientation bias and contamination, then filter Mutect2 calls. Conda env: wes.
set -euo pipefail
REF=ref/Homo_sapiens_assembly38.fasta
COMMON=ref/resources/small_exac_common_3.hg38.vcf.gz

gatk LearnReadOrientationModel -I vcf/HCC1395.f1r2.tar.gz -O vcf/HCC1395.read-orientation-model.tar.gz

# Allele counts at common SNP sites: used for contamination here and re-used for BAF in step 11
for s in HCC1395 HCC1395BL; do
  gatk GetPileupSummaries -I bam/$s.dedup.bam -V $COMMON -L $COMMON -O qc/$s.pileups.table
done
gatk CalculateContamination -I qc/HCC1395.pileups.table -matched qc/HCC1395BL.pileups.table \
  --tumor-segmentation qc/HCC1395.segments.table -O qc/HCC1395.contamination.table

gatk FilterMutectCalls -R $REF -V vcf/HCC1395.unfiltered.vcf.gz \
  --ob-priors vcf/HCC1395.read-orientation-model.tar.gz \
  --contamination-table qc/HCC1395.contamination.table \
  --tumor-segmentation qc/HCC1395.segments.table \
  -O vcf/HCC1395.filtered.vcf.gz

bcftools view -f PASS vcf/HCC1395.filtered.vcf.gz -Oz -o vcf/HCC1395.pass.vcf.gz
bcftools index -t vcf/HCC1395.pass.vcf.gz
bcftools view -H vcf/HCC1395.pass.vcf.gz | wc -l   # 1,750 PASS variants
