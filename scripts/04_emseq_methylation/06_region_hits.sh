#!/usr/bin/env bash
# Every cytosine in every W22 gene interval and control region, per genotype and tissue.
#
# The region file holds the 16 W22 gene sets and the shuffled control regions (named ctrl_1 ..
# ctrl_6181 in the control file's order). One allc stream per sample, all regions in one pass.
# Gene and region methylation is computed from these hits downstream (06_per_gene_methylation,
# Supplemental Fig. S4).
#
# Input   $GI_WORK/04_emseq_methylation/allc/allc_<sample>_merged.tsv.gz   (04_merge_allc.sh)
#         $GI_WORK/02_gene_sets/w22/*.bed                                   (02_gene_sets)
#         $GI_WORK/04_emseq_methylation/metagene/shuffled_control.bed       (05_shuffled_control.sh)
# Output  $GI_WORK/04_emseq_methylation/violin/regions_with_control.bed
#         $GI_WORK/04_emseq_methylation/violin/<sample>_hits.tsv
#           set  region_id  context  mc  total   (one row per cytosine and region)
#
# Chromosome names: the allc tables use 1..10, the BED files chr1..chr10, so "chr" is removed.
set -euo pipefail
source "$(dirname "$0")/../../config.sh"
gi_env gi_reads

GS=$GI_WORK/02_gene_sets/w22
CTRL=$GI_WORK/04_emseq_methylation/metagene/shuffled_control.bed
ALLC=$GI_WORK/04_emseq_methylation/allc
WORK=$GI_WORK/04_emseq_methylation/violin
mkdir -p "$WORK"

REG=$WORK/regions_with_control.bed
: > "$REG.tmp"
for f in "$GS"/*.bed; do
    s=$(basename "$f" .bed)
    awk -v s="$s" 'BEGIN{OFS="\t"} {sub(/^chr/,"",$1); print $1,$2,$3,s,$4}' "$f" >> "$REG.tmp"
done
awk 'BEGIN{OFS="\t"} {sub(/^chr/,"",$1); print $1,$2,$3,"shuffled_control","ctrl_"NR}' \
    "$CTRL" >> "$REG.tmp"
sort -k1,1 -k2,2n "$REG.tmp" > "$REG"; rm -f "$REG.tmp"
echo "regions: $(wc -l < "$REG") rows, $(awk '$4=="shuffled_control"' "$REG" | wc -l) control"

for g in wt_15dap embryo mdr1_15dap wt_mature mdr1_mature dng102_mature leaf tassel; do
    hits="$WORK/${g}_hits.tsv"
    [ -s "$hits" ] && { echo "### skip $g (done)"; continue; }
    echo "### $g $(date '+%H:%M')"
    zcat "$ALLC/allc_${g}_merged.tsv.gz" \
        | awk 'BEGIN{OFS="\t"} $1!="unmapped"{print $1,$2-1,$2,$4,$5,$6}' \
        | bedtools intersect -a stdin -b "$REG" -wa -wb \
        | awk 'BEGIN{OFS="\t"} {print $10,$11,$4,$5,$6}' > "$hits.part"
    mv -f "$hits.part" "$hits"
    echo "    $(wc -l < "$hits") cytosine-region rows"
done
