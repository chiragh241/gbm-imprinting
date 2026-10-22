#!/usr/bin/env bash
# B73 AGPv4 gene coordinates, to match the coordinate space of the GSE120304 ChIP BEDs.
#
# Input   $GI_DATA/references/B73v4.gff3.gz   (00_references; MaizeGDB Zm00001d.2)
# Output  $GI_WORK/09_chip_leaf_ear/B73v4_genes.bed   chrom start end gene
set -u
source "$(dirname "$0")/../../config.sh"
D=$GI_WORK/09_chip_leaf_ear
mkdir -p "$D"
GFF=$GI_DATA/references/B73v4.gff3.gz
[ -s "$GFF" ] || { echo "!!! no v4 GFF"; exit 1; }

echo
echo "=== gene rows, and id format ==="
zcat "$GFF" | awk -F'\t' '$3=="gene"' | head -2
echo "  total gene rows: $(zcat "$GFF" | awk -F'\t' '$3=="gene"' | wc -l)"

# v4 gene BED: chrom start end gene
BED=$D/B73v4_genes.bed
zcat "$GFF" | awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene"{
    id=""
    n=split($9,a,";")
    for(i=1;i<=n;i++){
      if (a[i] ~ /^ID=/)      { id=a[i]; sub(/^ID=/,"",id); sub(/^gene:/,"",id) }
      else if (a[i] ~ /^gene_id=/ && id=="") { id=a[i]; sub(/^gene_id=/,"",id) }
    }
    if (id!="") print $1,$4-1,$5,id
  }' | sort -k1,1 -k2,2n > "$BED"
echo
echo "=== $BED ==="
head -3 "$BED"
echo "  rows: $(wc -l < "$BED")"
echo "  chroms: $(cut -f1 "$BED" | sort -u | tr '\n' ' ')"
echo "  Zm00001d ids: $(grep -c 'Zm00001d' "$BED")"
