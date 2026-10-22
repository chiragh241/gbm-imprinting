#!/usr/bin/env bash
# Fixed 100 bp metagene bins for Figure 2a-d and Supplemental Fig. S3.
#
#   bins  0-19  2 kb upstream flank, 20 x 100 bp
#   bins 20-39  2 kb inward from the TSS, 20 x 100 bp
#   bins 40-59  2 kb inward from the 3' end, 20 x 100 bp
#   bins 60-79  2 kb downstream flank, 20 x 100 bp
#
# A gene contributes to an inward bin only where the bin lies wholly inside the gene and does
# not cross the gene midpoint, so no base pair is counted twice and short genes drop out of the
# deeper bins; genes shorter than 200 bp are skipped. Bins 39 and 40 are not adjacent in the
# genome: the figure draws a break there. Bins are ordered 5' to 3' so both strands stack.
#
# Input   $GI_WORK/02_gene_sets/w22/{PEG,biparental}__{gbM,UM}.bed
#         $GI_WORK/04_emseq_methylation/metagene/shuffled_control.bed
# Output  $GI_WORK/04_emseq_methylation/metagene/metagene_bins_fixed.bed
#           chrom start end set gene bin
#         $GI_WORK/04_emseq_methylation/metagene/metagene_fixed_bin_counts.tsv
#           regions reaching each bin, per set
set -euo pipefail
source "$(dirname "$0")/../../config.sh"

GS=$GI_WORK/02_gene_sets/w22
WORK=$GI_WORK/04_emseq_methylation/metagene
SHUF=$WORK/shuffled_control.bed
mkdir -p "$WORK"

NUP=20; NIN=20; NDOWN=20; BINW=100
SETS="PEG__gbM biparental__gbM PEG__UM biparental__UM"

binify() {  # bedfile setname
  awk -v S="$2" -v NUP=$NUP -v NIN=$NIN -v ND=$NDOWN -v BW=$BINW '
    BEGIN{OFS="\t"}
    {
      chrom=$1; sub(/^chr/,"",chrom); s=$2; e=$3; name=$4; strand=$6
      L=e-s; if (L < 2*BW) next            # too short for even one inward bin each side
      mid=s+int(L/2)
      n=NUP+2*NIN+ND
      for (i=0;i<NUP;i++){ a=s-(NUP-i)*BW; b=s-(NUP-i-1)*BW; if(a<0) continue
                           k=i; if(strand=="-") k=n-1-k; print chrom,a,b,S,name,k }
      for (i=0;i<NIN;i++){ a=s+i*BW; b=s+(i+1)*BW; if (b>mid) break
                           k=NUP+i; if(strand=="-") k=n-1-k; print chrom,a,b,S,name,k }
      for (i=0;i<NIN;i++){ a=e-(NIN-i)*BW; b=e-(NIN-i-1)*BW; if (a<mid) continue
                           k=NUP+NIN+i; if(strand=="-") k=n-1-k; print chrom,a,b,S,name,k }
      for (i=0;i<ND;i++) { a=e+i*BW; b=e+(i+1)*BW
                           k=NUP+2*NIN+i; if(strand=="-") k=n-1-k; print chrom,a,b,S,name,k }
    }' "$1"
}

BINS=$WORK/metagene_bins_fixed.bed
: > "$BINS.tmp"
for s in $SETS; do
    binify "$GS/$s.bed" "$s" >> "$BINS.tmp"
    echo "binned $s ($(wc -l < "$GS/$s.bed") genes)"
done
binify "$SHUF" "shuffled_control" >> "$BINS.tmp"
sort -k1,1 -k2,2n "$BINS.tmp" > "$BINS"; rm -f "$BINS.tmp"
echo "bins: $(wc -l < "$BINS") rows"

# The header is added after the sort, which would otherwise file it among the data rows.
awk 'BEGIN{OFS="\t"} {n[$4"\t"$6]++} END{for(k in n) print k,n[k]}' "$BINS" \
    | sort -k1,1 -k2,2n \
    | { printf 'set\tbin\tn_regions\n'; cat; } > "$WORK/metagene_fixed_bin_counts.tsv"
