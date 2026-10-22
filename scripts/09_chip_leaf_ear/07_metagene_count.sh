#!/usr/bin/env bash
# Count Ricci 2019 (GSE120304) H3K27me3 and input reads into the metagene bins.
#
# Same streaming pattern as the Figure 1E counting step: reads are -a, the bins are -b, so
# bedtools holds only the bin side. Libraries run one at a time; 1.2M bins as -b is a large
# interval tree and two of those at once is how the VM gets wedged.
#
# Library totals are not recounted: 04_count.sh already wrote them (counts/*.total.txt).
# Every bin a read overlaps is counted, so a read across a bin edge counts in both bins.
#
# Input   $GI_DATA/ricci2019/<lib>.bed.gz, $GI_WORK/09_chip_leaf_ear/metagene_bins_v4.bed
# Output  $GI_WORK/09_chip_leaf_ear/bincounts/<lib>.bincounts.tsv   set gene scheme bin reads
set -u
source "$(dirname "$0")/../../config.sh"
gi_env gi_analysis
D=$GI_WORK/09_chip_leaf_ear
RAW=$GI_DATA/ricci2019
TOT=$D/counts
BINS=$D/metagene_bins_v4.bed
OUT=$D/bincounts
mkdir -p "$OUT"

[ -s "$BINS" ] || { echo "!!! $BINS missing, run 06_metagene_bins.py first"; exit 1; }
echo "=== start $(date)"
echo "bins: $(wc -l < "$BINS") rows"

for lib in leaf_H3K27me3_rep1 leaf_H3K27me3_rep2 leaf_input_rep1 leaf_input_rep2 \
           ear_H3K27me3_rep1 ear_H3K27me3_rep2 ear_input_rep1 ; do
  f=$RAW/${lib}.bed.gz
  o=$OUT/${lib}.bincounts.tsv
  [ -s "$f" ] || { echo "### $lib: no file, skip"; continue; }
  [ -s "$o" ] && { echo "### skip $lib (done)"; continue; }
  [ -s "$TOT/${lib}.total.txt" ] || { echo "!!! $lib has no total.txt"; exit 1; }
  echo "### $lib  $(date '+%H:%M')  total reads $(cat "$TOT/${lib}.total.txt")"
  # bins cols after the 6 read cols: 7 chrom 8 start 9 end 10 set 11 gene 12 scheme 13 bin
  zcat "$f" \
    | nice -n 19 bedtools intersect -a stdin -b "$BINS" -wa -wb \
    | awk 'BEGIN{OFS="\t"} {n[$10"\t"$11"\t"$12"\t"$13]++}
           END{for(k in n) print k,n[k]}' \
    > "$o.part"
  mv -f "$o.part" "$o"
  echo "    cells: $(wc -l < "$o")"
done

echo "=== done $(date)"
ls -la "$OUT"
