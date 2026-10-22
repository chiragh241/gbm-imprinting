#!/usr/bin/env bash
# H3K27me3 ChIP-seq and input of B73 seedling leaf and developing ear (Supplemental Fig. S1).
#
# Source: Ricci et al. 2019 Nature Plants (Schmitz lab), GEO GSE120304, supplementary files
# <GSM>_*.bed.gz (mapped reads, B73 AGPv4):
#   GSM3397987 leaf H3K27me3 rep1   GSM3397988 leaf H3K27me3 rep2
#   GSM3398001 leaf input rep1      GSM3398002 leaf input rep2
#   GSM3398003 ear H3K27me3 rep1    GSM3398004 ear H3K27me3 rep2   GSM3398012 ear input rep1
#
# Files are MAPPED READS in BED (AGPv4), not peak calls, so enrichment is computed here rather
# than taken from someone's threshold. Input samples are downloaded too: without input
# normalisation a mappability or copy-number difference between gbM and UM genes could look
# like an H3K27me3 difference, which is exactly the comparison this figure makes.
#
# Output  $GI_DATA/ricci2019/<label>.bed.gz
set -u
source "$(dirname "$0")/../../config.sh"
RAW=$GI_DATA/ricci2019
mkdir -p "$RAW"

declare -A S=(
  [GSM3397987]="leaf_H3K27me3_rep1"
  [GSM3397988]="leaf_H3K27me3_rep2"
  [GSM3398001]="leaf_input_rep1"
  [GSM3398002]="leaf_input_rep2"
  [GSM3398003]="ear_H3K27me3_rep1"
  [GSM3398004]="ear_H3K27me3_rep2"
  [GSM3398012]="ear_input_rep1"
)

base=https://ftp.ncbi.nlm.nih.gov/geo/samples
bad=0
for g in "${!S[@]}"; do
  label=${S[$g]}
  out=$RAW/${label}.bed.gz
  if [ -s "$out" ]; then echo "### skip $label (have it)"; continue; fi
  pre=${g%???}nnn
  fn=$(curl -s --max-time 120 "$base/$pre/$g/suppl/" \
        | grep -oE "${g}_[^\"<>]+\.bed\.gz" | sort -u | head -1)
  if [ -z "$fn" ]; then echo "!!! $g: no bed.gz found"; bad=1; continue; fi
  echo "### $label <- $fn $(date '+%H:%M')"
  curl -sf --max-time 3600 -o "$out.part" "$base/$pre/$g/suppl/$fn" \
    && gzip -t < "$out.part" \
    && mv -f "$out.part" "$out" \
    && echo "    $(du -h "$out" | cut -f1)" \
    || { echo "!!! download failed $g"; rm -f "$out.part"; bad=1; }
done

echo
echo "=== downloaded ==="
ls -la "$RAW" | awk '{printf "  %-34s %s\n",$9,$5}'
echo
echo "=== sanity: first lines of one file ==="
zcat "$RAW/leaf_H3K27me3_rep1.bed.gz" 2>/dev/null | head -3
echo "=== chrom names present ==="
zcat "$RAW/leaf_H3K27me3_rep1.bed.gz" 2>/dev/null | cut -f1 | head -200000 | sort -u | head -15
# a missing library would otherwise be skipped by the counting steps and change the pooled values
[ "$bad" = 0 ] || { echo "!!! not every library downloaded; run again"; exit 1; }
