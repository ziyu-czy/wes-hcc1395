#!/bin/bash
# Step 4: somatic SNV/indel calling with Mutect2 (tumor vs matched normal). Conda env: wes (gatk4).
# On Apple Silicon this runs single-threaded (~5 h for the whole exome); caffeinate keeps the Mac awake.
set -euo pipefail
mkdir -p vcf
caffeinate -i gatk --java-options "-Xmx16g" Mutect2 \
  -R ref/Homo_sapiens_assembly38.fasta \
  -I bam/HCC1395.dedup.bam -I bam/HCC1395BL.dedup.bam -normal HCC1395BL \
  --germline-resource ref/resources/af-only-gnomad.hg38.vcf.gz \
  --panel-of-normals ref/resources/1000g_pon.hg38.vcf.gz \
  -L ref/exome_targets.hg38.bed \
  --f1r2-tar-gz vcf/HCC1395.f1r2.tar.gz \
  -O vcf/HCC1395.unfiltered.vcf.gz
