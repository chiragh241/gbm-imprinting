#!/usr/bin/env bash
# Nearest DMR to each gene, with its distance relative to the gene's orientation.
#
#   bedtools closest -D a -t first   (bedtools 2.31.1)
#     -D a     negative when the DMR is upstream, positive downstream, 0 when it overlaps
#     -t first one DMR when several are equally close
#
# Input   $GI_WORK/07_dmr_distances/genes.sorted.bed, dmrs.sorted.bed   (01_prep_beds.py)
# Output  $GI_WORK/07_dmr_distances/closest.tsv
set -euo pipefail
source "$(dirname "$0")/../../config.sh"
gi_env gi_analysis

W=$GI_WORK/07_dmr_distances
bedtools closest -a "$W/genes.sorted.bed" -b "$W/dmrs.sorted.bed" -D a -t first > "$W/closest.tsv"
echo "closest.tsv: $(wc -l < "$W/closest.tsv") genes"
