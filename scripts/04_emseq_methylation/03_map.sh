#!/usr/bin/env bash
# Map each EM-seq library to W22 v2 and call methylation with the methylpy paired-end pipeline.
#
#   methylpy paired-end-pipeline --num-procs 8 --remove-clonal True --sort-mem 2G
#     cutadapt 5.2 with methylpy's default adapters (read 1 AGATCGGAAGAGCACACGTCTGAAC,
#     read 2 AGATCGGAAGAGCGTCGTGTAGGGA) and minimum read length 30; bowtie2 2.5.5;
#     clonal reads removed with Picard 3.4.0 MarkDuplicates; MAPQ >= 30; per-cytosine
#     methylated and total read counts written as an allc table.
#   methylpy is patched to run with these cutadapt and Picard versions (envs/patches/methylpy.patch,
#   applied by envs/setup_envs.sh).
#
# Each replicate is mapped and deduplicated on its own, and replicates are summed per cytosine
# afterwards (04_merge_allc.sh): deduplicating after merging BAMs would treat reads that collide
# across biological replicates as PCR duplicates.
#
# Input   $GI_DATA/reads/emseq_w22/<tag>_{1,2}.fastq.gz   (01_download_reads.sh)
#         $GI_WORK/04_emseq_methylation/genome/W22v2_{f,r}  (02_bisulfite_reference.sh)
# Output  $GI_WORK/04_emseq_methylation/allc/allc_<tag>.tsv.gz, one per library
#
# About 1 hour of wall time per Gbase of reads on 8 cores: this step took several weeks here.
# The merged methylation calls it leads to are deposited on Zenodo (see README), so the
# downstream steps can start from there.
set -u
source "$(dirname "$0")/../../config.sh"
gi_env gi_reads

FQ=$GI_DATA/reads/emseq_w22
W=$GI_WORK/04_emseq_methylation
G=$W/genome
OUT=$W/allc
MAPW=$W/mapwork
# Temporary files go to the work disk: Picard's sort spill does not fit a RAM-backed /tmp.
# TMPDIR covers sort; java.io.tmpdir covers Picard, which ignores TMPDIR.
TMP=$W/tmp
mkdir -p "$OUT" "$MAPW" "$TMP"
export TMPDIR="$TMP"
export _JAVA_OPTIONS="-Djava.io.tmpdir=$TMP"

FA=$GI_DATA/references/W22v2.fa
PDIR=$(dirname "$(find "$CONDA_PREFIX/share" -name 'picard.jar' 2>/dev/null | head -1)")
[ -n "$PDIR" ] || { echo "!!! picard.jar not found"; exit 1; }

FAILED=0
while IFS=$'\t' read -r run tag bytes; do
    [ "$run" = "run" ] && continue
    if [ -s "$OUT/allc_${tag}.tsv.gz" ]; then echo "### skip $tag (allc exists)"; continue; fi
    r1="$FQ/${tag}_1.fastq.gz"; r2="$FQ/${tag}_2.fastq.gz"
    [ -s "$r1" ] && [ -s "$r2" ] || { echo "!!! missing fastq for $tag"; FAILED=1; continue; }
    echo "### $tag  $(date)"
    ( cd "$MAPW" && methylpy paired-end-pipeline \
        --read1-files "$r1" --read2-files "$r2" --sample "$tag" \
        --forward-ref "$G/W22v2_f" --reverse-ref "$G/W22v2_r" --ref-fasta "$FA" \
        --num-procs 8 --remove-clonal True --path-to-output "$OUT" \
        --path-to-picard "$PDIR" --java-options=-Xmx8g \
        --sort-mem 2G ) || { echo "!!! methylpy failed on $tag"; FAILED=1; }
done < "$GI_SHIPPED/emseq/runs.tsv"
exit $FAILED
