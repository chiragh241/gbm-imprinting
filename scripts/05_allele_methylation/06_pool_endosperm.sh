#!/usr/bin/env bash
# The six deduplicated endosperm libraries merged into one BAM, for SNV calling only.
#
# Variants are called once on all six libraries pooled, both cross directions together, so the
# variant set is direction-neutral; each library is then phased from its own BAM (10_asm.sh).
# The libraries are merged at full depth (51 to 81 million mapped reads each, a 1.58-fold
# spread): all six share the same two parents, so depth differences cannot change which base is
# called B73 or W22, and downsampling to the shallowest library would push sites under the 15x
# coverage floor of 08_snv.sh. Bisulfite conversion is uniform across the six (mCHH 1.72 to
# 1.94%), so no library is over-weighted on that account either.
#
# Input   $GI_WORK/05_allele_methylation/dedup/<library>.dedup.bam   (05_sort_dedup.sh)
# Output  $GI_WORK/05_allele_methylation/combined/endosperm_combined.bam (+ .bai)
set -euo pipefail
source "$(dirname "$0")/../../config.sh"
gi_env gi_bsseeker2

W=$GI_WORK/05_allele_methylation
OUT=$W/combined
mkdir -p "$OUT"
BAM=$OUT/endosperm_combined.bam
[ -f "$OUT/.combine_done" ] && { echo "already combined"; exit 0; }

inputs=()
for l in $(awk -F'\t' '$2=="endosperm"{print $1}' "$W/libs.tsv"); do
    inputs+=("$W/dedup/${l}.dedup.bam")
    echo "$l: $(samtools idxstats "$W/dedup/${l}.dedup.bam" | awk '{s+=$3} END{print s}') mapped reads"
done
samtools merge -@ "$GI_THREADS" -f "$BAM" "${inputs[@]}"
samtools index -@ "$GI_THREADS" "$BAM"
samtools quickcheck "$BAM" && touch "$OUT/.combine_done"
