#!/usr/bin/env bash
# Coordinate sort and duplicate removal of each library (samtools 1.18):
#   collate -> fixmate -m -> sort -> markdup -r -> index
# markdup needs the mate-score tags fixmate adds, which needs name-collated input.
#
# Input   $GI_WORK/05_allele_methylation/map/<library>.bam       (04_map.sh)
# Output  $GI_WORK/05_allele_methylation/dedup/<library>.dedup.bam (+ .bai, markdup stats)
set -u
source "$(dirname "$0")/../../config.sh"
gi_env gi_bsseeker2

W=$GI_WORK/05_allele_methylation
OUT=$W/dedup
TMP=$W/tmp
T=${GI_THREADS}
mkdir -p "$OUT" "$TMP"

FAILED=0
while IFS=$'\t' read -r lib tissue direction r1 r2; do
    [ -f "$OUT/${lib}.dedup_done" ] && { echo "### skip $lib"; continue; }
    in=$W/map/${lib}.bam
    col=$TMP/${lib}.collate.bam; fix=$TMP/${lib}.fixmate.bam; crd=$TMP/${lib}.coord.bam
    ded=$OUT/${lib}.dedup.bam
    echo "### $lib $(date)"
    samtools collate -@ "$T" -o "$col" "$in" "$TMP/${lib}.col" \
      && samtools fixmate -@ "$T" -m "$col" "$fix" && rm -f "$col" \
      && samtools sort -@ "$T" -m 1G -T "$TMP/${lib}.srt" -o "$crd" "$fix" && rm -f "$fix" \
      && samtools markdup -@ "$T" -r -T "$TMP/${lib}.mkd" -f "$OUT/${lib}.markdup_stats.txt" \
             "$crd" "$ded" && rm -f "$crd" \
      && samtools index -@ "$T" "$ded" && samtools quickcheck "$ded" \
      && touch "$OUT/${lib}.dedup_done" \
      || { echo "!!! $lib failed"; FAILED=1; }
done < "$W/libs.tsv"
exit $FAILED
