#!/usr/bin/env bash
# Forward and reverse bisulfite-converted W22 v2 references for methylpy, and the W22 v2 index.
#
# Input   $GI_DATA/references/W22v2.fa   (00_references/01_download_references.sh)
# Output  $GI_DATA/references/W22v2.fa.fai
#         $GI_WORK/04_emseq_methylation/genome/W22v2_f.*, W22v2_r.*   (bowtie2 indexes)
set -euo pipefail
source "$(dirname "$0")/../../config.sh"
gi_env gi_reads

FA=$GI_DATA/references/W22v2.fa
G=$GI_WORK/04_emseq_methylation/genome
mkdir -p "$G"

[ -s "$FA.fai" ] || samtools faidx "$FA"
if [ ! -f "$G/.build_done" ]; then
    methylpy build-reference \
        --input-files "$FA" \
        --output-prefix "$G/W22v2" \
        --aligner bowtie2 \
        --num-procs "$GI_THREADS"
    touch "$G/.build_done"
fi
ls -la "$G"
