#!/usr/bin/env bash
# CG methylated and total reads per gene and metagene bin, in the four tissues of Figure 2a-d.
#
# Methylation of a bin is summed over its CG cytosines (second base of the trinucleotide G).
# Keeping the gene lets Figure 2a-d pool over genes and report how many genes reach each bin.
#
# Input   $GI_WORK/04_emseq_methylation/allc/allc_<tissue>_merged.tsv.gz
#         $GI_WORK/04_emseq_methylation/metagene/metagene_bins_fixed.bed   (07_metagene_bins.sh)
# Output  $GI_WORK/04_emseq_methylation/metagene/<tissue>_metagene_cg_fixed_pergene.tsv
#           set  gene  bin  mc  total        tissues: wt_15dap, embryo, leaf, tassel
set -euo pipefail
source "$(dirname "$0")/../../config.sh"
gi_env gi_reads

ALLC=$GI_WORK/04_emseq_methylation/allc
WORK=$GI_WORK/04_emseq_methylation/metagene
BINS=$WORK/metagene_bins_fixed.bed

for t in wt_15dap embryo leaf tassel; do
    out="$WORK/${t}_metagene_cg_fixed_pergene.tsv"
    [ -s "$out" ] && { echo "### skip $t (done)"; continue; }
    echo "### $t $(date '+%H:%M')"
    # bins columns after the cytosine's six: 7 chrom 8 start 9 end 10 set 11 gene 12 bin
    zcat "$ALLC/allc_${t}_merged.tsv.gz" \
        | awk 'BEGIN{OFS="\t"} $1!="unmapped" && substr($4,2,1)=="G" {print $1,$2-1,$2,$4,$5,$6}' \
        | bedtools intersect -a stdin -b "$BINS" -wa -wb \
        | awk 'BEGIN{OFS="\t"} {k=$10"\t"$11"\t"$12; mc[k]+=$5; tot[k]+=$6}
               END{for(k in mc) print k,mc[k],tot[k]}' > "$out.body"
    { printf 'set\tgene\tbin\tmc\ttotal\n'; sort -k1,1 -k2,2 -k3,3n "$out.body"; } > "$out.part"
    mv -f "$out.part" "$out"
    rm -f "$out.body"
    echo "    $(( $(wc -l < "$out") - 1 )) set/gene/bin rows"
done
