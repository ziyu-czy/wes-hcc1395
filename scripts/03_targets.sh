#!/bin/bash
# Step 3: build the exome target BED (Agilent SureSelect V6+UTR, S07604624, hg38 liftover from the SEQC2 FTP).
# Keep only primary chromosomes, sort in reference order, merge overlaps. Conda env: wes (bedtools).
set -euo pipefail
cd ref
zcat < S07604624_Covered_human_all_v6_plus_UTR.liftover.to.hg38.bed6.gz \
  | awk '$1 ~ /^chr([0-9]+|X|Y)$/' | cut -f1-3 \
  | bedtools sort -faidx Homo_sapiens_assembly38.fasta.fai \
  | bedtools merge > exome_targets.hg38.bed
wc -l exome_targets.hg38.bed   # expected: 294,866 regions, ~90.5 Mb
