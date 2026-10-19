#!/bin/bash
# Step 10: copy-number ratios and segmentation with CNVkit (tumor vs matched normal).
# Conda env with cnvkit (and R DNAcopy). Run from the project root.
set -euo pipefail
mkdir -p cnv results/figures

# Accessible (non-N) regions of the reference, primary chromosomes only
cnvkit.py access ref/Homo_sapiens_assembly38.fasta -o cnv/access.all.bed
grep -E '^chr([0-9]+|X|Y)\b' cnv/access.all.bed > cnv/access.hg38.bed

caffeinate -i cnvkit.py batch bam/HCC1395.dedup.bam \
  --normal bam/HCC1395BL.dedup.bam \
  --targets ref/exome_targets.hg38.bed \
  --access cnv/access.hg38.bed \
  --fasta ref/Homo_sapiens_assembly38.fasta \
  --output-reference cnv/HCC1395BL.reference.cnn \
  --output-dir cnv/ \
  --drop-low-coverage \
  --scatter --diagram \
  -p 8

cnvkit.py metrics cnv/HCC1395.dedup.cnr -s cnv/HCC1395.dedup.cns   # mad 0.33, bivar 0.35 (noisy)
cp cnv/HCC1395.dedup-scatter.png results/figures/cnvkit_scatter.png
# Note: HCC1395.dedup.call.cns assumes a diploid, pure sample -> not used (tumor is near-triploid).
