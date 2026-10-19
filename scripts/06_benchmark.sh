#!/bin/bash
# Step 6: benchmark PASS calls against the SEQC2 high-confidence somatic truth set (v1.2.1).
# Evaluated only inside exome targets ∩ SEQC2 high-confidence regions. Conda env: wes.
# isec output: 0000 = only ours (FP), 0001 = only truth (FN), 0002/0003 = shared (TP).
set -euo pipefail
REF=ref/Homo_sapiens_assembly38.fasta
mkdir -p bench results/tables

bedtools intersect -a ref/exome_targets.hg38.bed -b truth/High-Confidence_Regions_v1.2.bed \
  | bedtools sort -faidx ref/Homo_sapiens_assembly38.fasta.fai | bedtools merge > truth/eval_regions.bed

printf "type\tTP\tFP\tFN\tprecision\trecall\n" > results/tables/benchmark.tsv
for t in snps indels; do
  if [ $t = snps ]; then TRUTH=truth/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz
  else TRUTH=truth/high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz; fi
  bcftools index -f -t $TRUTH
  bcftools view -f PASS -v $t -R truth/eval_regions.bed vcf/HCC1395.filtered.vcf.gz \
    | bcftools norm -f $REF -m -any -Oz -o bench/ours_$t.vcf.gz
  bcftools index -f -t bench/ours_$t.vcf.gz
  bcftools view -v $t -R truth/eval_regions.bed $TRUTH \
    | bcftools norm -f $REF -m -any -Oz -o bench/truth_$t.vcf.gz
  bcftools index -f -t bench/truth_$t.vcf.gz
  bcftools isec -c none -p bench/isec_$t bench/ours_$t.vcf.gz bench/truth_$t.vcf.gz
  FP=$(grep -vc '^#' bench/isec_$t/0000.vcf || true)
  FN=$(grep -vc '^#' bench/isec_$t/0001.vcf || true)
  TP=$(grep -vc '^#' bench/isec_$t/0002.vcf || true)
  awk -v t=$t -v tp=$TP -v fp=$FP -v fn=$FN \
    'BEGIN{printf "%s\t%d\t%d\t%d\t%.3f\t%.3f\n", t, tp, fp, fn, tp/(tp+fp), tp/(tp+fn)}' \
    | tee -a results/tables/benchmark.tsv
done
