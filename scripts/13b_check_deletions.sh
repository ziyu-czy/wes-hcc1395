#!/bin/bash
# Step 13b: confirm homozygous deletions directly from read depth, exon by exon (MANE Select exons, GENCODE 50).
# A VCF cannot show this: where there are no reads there is nothing to call. Conda env: wes.
# Result: CDKN2A - all exons ~0x in tumor (normal 47-284x) -> whole-gene homozygous deletion
#         PTEN   - 5' exons retained, 3' part (from ~chr10:87.957 Mb to the gene end) ~0x -> partial homozygous deletion
set -euo pipefail
GTF=ref/gencode.v50.annotation.gtf.gz

for g in CDKN2A PTEN; do
  # GTF is 1-based, BED is 0-based -> start - 1
  zcat < $GTF | grep "gene_name \"$g\";" | grep 'MANE_Select' \
  | awk -F'\t' -v g=$g 'BEGIN{OFS="\t"} $3=="exon" {match($9, /exon_number [0-9]+/);
        print $1, $4-1, $5, g "_exon" substr($9, RSTART+12, RLENGTH-12)}' > cnv/${g}_exons.bed

  samtools bedcov cnv/${g}_exons.bed bam/HCC1395.dedup.bam bam/HCC1395BL.dedup.bam \
  | awk 'BEGIN{OFS="\t"; print "chrom","start","end","exon","tumor_DP","normal_DP"}
         {L=$3-$2; printf "%s\t%d\t%d\t%s\t%.1f\t%.1f\n", $1, $2, $3, $4, $5/L, $6/L}' \
  | tee results/tables/${g}_exon_depth.tsv | column -t
done
