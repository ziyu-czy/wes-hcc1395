#!/bin/bash
# Step 1: QC and alignment of tumor (HCC1395) and normal (HCC1395BL) WES reads.
# Run from the project root. Conda env: wes (fastqc, bwa, samtools).
set -euo pipefail
REF=ref/Homo_sapiens_assembly38.fasta   # GATK hg38 bundle; bwa index + .alt file in ref/

mkdir -p qc bam
fastqc --outdir qc fastq/SRR7890883_1.fastq.gz fastq/SRR7890883_2.fastq.gz \
                   fastq/SRR7890874_1.fastq.gz fastq/SRR7890874_2.fastq.gz

align() {
  local sample=$1 run=$2 rgid=$3 lib=$4
  # -K fixes the batch size so results are reproducible regardless of thread count
  # -Y soft-clips supplementary alignments; the read group SM tag names the sample for Mutect2
  bwa mem -t 8 -K 100000000 -Y \
    -R "@RG\tID:${rgid}\tSM:${sample}\tLB:${lib}\tPL:ILLUMINA\tPU:${run}" \
    $REF fastq/${run}_1.fastq.gz fastq/${run}_2.fastq.gz \
    2> bam/${sample}.bwa.log \
  | samtools sort -@ 4 -m 2G -o bam/${sample}.sorted.bam -
  samtools index bam/${sample}.sorted.bam
  samtools flagstat bam/${sample}.sorted.bam > bam/${sample}.flagstat.txt
}

align HCC1395   SRR7890883 HCC1395_T1   WES_IL_T_1
align HCC1395BL SRR7890874 HCC1395BL_N1 WES_IL_N_1
echo "ALL DONE"
