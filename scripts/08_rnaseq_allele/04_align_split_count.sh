#!/usr/bin/env bash
# Per sample: HISAT2 to the N-masked B73 v5, SNPsplit into B73 (genome1) and Mo17 (genome2)
# alleles, and featureCounts of fragments per gene over exons. Reads are not trimmed.
#
#   hisat2 --no-mixed --no-discordant --no-unal          (HISAT2 2.2.1)
#   SNPsplit --paired                                     (SNPsplit 0.6.0)
#   featureCounts -p --countReadPairs -s 0 -F SAF         (Subread 2.0.8)
#
# Input   $GI_DATA/reads/rnaseq_xu2022/<run>_{1,2}.fastq.gz   (01_download_reads.sh)
#         data/rnaseq/samplesheet.tsv
#         $GI_WORK/08_rnaseq_allele/{index,snps,annot}          (02, 03)
# Output  $GI_WORK/08_rnaseq_allele/counts/<run>.B73.txt, <run>.Mo17.txt
#
# Files are named by SRA run, not by sample: the sample names differ only in case (BM_r1 is
# B73 x Mo17, bM_r1 the zmros1ab cross), so on a case-insensitive file system they would be one
# file. Once a sample is counted its alignments are deleted (about 23 GB a sample); the
# per-allele BAMs, the counts and the logs are kept.
#
# HISAT2's spliced alignment is not deterministic between runs: rerun, these counts differ from
# the paper's by one fragment in about 1% of genes. The paper's count tables are on Zenodo (see
# "Skipping the read mapping" in the README).
set -euo pipefail
source "$(dirname "$0")/../../config.sh"
gi_env gi_reads

W=$GI_WORK/08_rnaseq_allele
IDX=$W/index/B73v5_Nmasked
SNP=$W/snps/snps.snpsplit.txt
SAF=$W/annot/exons.saf
FQ=$GI_DATA/reads/rnaseq_xu2022
BAM=$W/bam; SPLIT=$W/split; CNT=$W/counts; LOG=$W/logs
mkdir -p "$BAM" "$SPLIT" "$CNT" "$LOG"

tail -n +2 "$GI_SHIPPED/rnaseq/samplesheet.tsv" | while IFS=$'\t' read -r sample cross mat pat geno rep run; do
    if [ -f "$CNT/${run}.done" ]; then
        echo "skip $sample ($run)"; continue
    fi
    echo ">>> $sample ($run) $(date)"
    B=$BAM/${run}.bam
    if [ ! -s "$B" ]; then
        hisat2 -p "$GI_THREADS" --no-mixed --no-discordant --no-unal -x "$IDX" \
            -1 "$FQ/${run}_1.fastq.gz" -2 "$FQ/${run}_2.fastq.gz" 2> "$LOG/hisat2_${run}.log" \
            | samtools view -bS -@ 2 - > "$B.part"
        mv -f "$B.part" "$B"
    fi
    SNPsplit --snp_file "$SNP" --paired -o "$SPLIT" "$B" > "$LOG/snpsplit_${run}.log" 2>&1
    featureCounts -p --countReadPairs -s 0 -T 8 -F SAF -a "$SAF" -o "$CNT/${run}.B73.txt" \
        "$SPLIT/${run}.genome1.bam" > "$LOG/fc_${run}_B73.log" 2>&1
    featureCounts -p --countReadPairs -s 0 -T 8 -F SAF -a "$SAF" -o "$CNT/${run}.Mo17.txt" \
        "$SPLIT/${run}.genome2.bam" > "$LOG/fc_${run}_Mo17.log" 2>&1
    rm -f "$B" "$BAM/${run}.sortedByName.bam" "$SPLIT/${run}.allele_flagged.bam" \
        "$SPLIT/${run}.unassigned.bam"
    touch "$CNT/${run}.done"
done
