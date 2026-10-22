#!/usr/bin/env bash
# Count ChIP and input reads over the gene bodies of Supplemental Fig. S1.
#
# The reads file is streamed as -a and the 6,055 gene intervals held as -b, so bedtools keeps
# only the small side in memory. That is the same pattern the methylome steps use and avoids
# both a 3.6 GB sort and an OOM.
#
# Total mapped reads per library are counted in the same pass via tee, so normalisation does
# not cost a second decompression of every file.
#
# Input   $GI_DATA/ricci2019/<lib>.bed.gz, $GI_WORK/09_chip_leaf_ear/fig1e_gene_intervals.bed
# Output  $GI_WORK/09_chip_leaf_ear/counts/<lib>.counts.tsv (gene, set, reads), <lib>.total.txt
set -u
source "$(dirname "$0")/../../config.sh"
gi_env gi_analysis
RAW=$GI_DATA/ricci2019
REG=$GI_WORK/09_chip_leaf_ear/fig1e_gene_intervals.bed
OUT=$GI_WORK/09_chip_leaf_ear/counts
mkdir -p "$OUT"

[ -s "$REG" ] || { echo "!!! $REG missing, run 03 first"; exit 1; }

for lib in leaf_H3K27me3_rep1 leaf_H3K27me3_rep2 leaf_input_rep1 leaf_input_rep2 \
           ear_H3K27me3_rep1 ear_H3K27me3_rep2 ear_input_rep1 ; do
  f=$RAW/${lib}.bed.gz
  o=$OUT/${lib}.counts.tsv
  t=$OUT/${lib}.total.txt
  [ -s "$f" ] || { echo "### $lib: no file, skip"; continue; }
  [ -s "$o" ] && [ -s "$t" ] && { echo "### skip $lib (counted)"; continue; }
  echo "### $lib $(date '+%H:%M')"
  zcat "$f" \
    | tee >(wc -l > "$t.part") \
    | nice -n 19 bedtools intersect -a stdin -b "$REG" -wa -wb \
    | awk 'BEGIN{OFS="\t"} {n[$10"\t"$11]++} END{for(k in n) print k,n[k]}' \
    > "$o.part"
  wait
  mv -f "$o.part" "$o"; mv -f "$t.part" "$t"
  echo "    genes hit: $(wc -l < "$o")   total reads: $(cat "$t")"
done

echo
echo "=== done $(date) ==="
ls -la "$OUT" | awk '{printf "  %-38s %s\n",$9,$5}'
