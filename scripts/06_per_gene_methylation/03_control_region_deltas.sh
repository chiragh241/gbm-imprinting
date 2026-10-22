#!/usr/bin/env bash
# CG methylation of each shuffled control region, per genotype, for the paired mutant deltas.
#
# Weighted over the region's CG cytosines (second base of the trinucleotide G), with the same
# filter as the genes: at least 5 covered CG sites and at least 10 reads in the region.
#
# Input   $GI_WORK/04_emseq_methylation/violin/<genotype>_hits.tsv   (04_emseq_methylation/06_region_hits.sh)
# Output  $GI_WORK/06_per_gene_methylation/ctrl_<genotype>.tsv   region_id  CG methylation
set -euo pipefail
source "$(dirname "$0")/../../config.sh"

V=$GI_WORK/04_emseq_methylation/violin
OUT=$GI_WORK/06_per_gene_methylation
mkdir -p "$OUT"

for g in wt_15dap mdr1_15dap wt_mature mdr1_mature dng102_mature; do
    out=$OUT/ctrl_${g}.tsv
    [ -s "$out" ] && { echo "### $g: done"; continue; }
    awk -F'\t' '
        $1 == "shuffled_control" && substr($3, 2, 1) == "G" { mc[$2] += $4; tot[$2] += $5; n[$2]++ }
        END { for (r in tot) if (n[r] >= 5 && tot[r] >= 10) printf "%s\t%.6f\n", r, mc[r]/tot[r] }
    ' "$V/${g}_hits.tsv" > "$out.part"
    mv -f "$out.part" "$out"
    echo "### $g: $(wc -l < "$out") control regions"
done
