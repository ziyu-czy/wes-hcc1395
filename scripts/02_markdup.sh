#!/bin/bash
# Step 2: mark PCR/optical duplicates. Conda env: wes (gatk4).
# BQSR was deliberately skipped (modern data; small effect on somatic SNV calls) - see README.
set -euo pipefail
mkdir -p qc
for s in HCC1395 HCC1395BL; do
  gatk --java-options "-Xmx16g" MarkDuplicates \
    -I bam/$s.sorted.bam -O bam/$s.dedup.bam -M qc/$s.dup_metrics.txt \
    --READ_NAME_REGEX null --CREATE_INDEX true
done
echo "ALL DONE"
