#!/usr/bin/env bash
# Shuffled control regions: the W22 gbM-biparental gene intervals placed at random over
# chromosomes 1 to 10, so the control matches that set in number and length but not in sequence.
#
#   bedtools shuffle -noOverlapping -seed 42, otherwise unconstrained
#
# Unconstrained on purpose: the control shows that maternal CG
# demethylation is not a general feature of the genome, so shuffled regions landing in
# methylated TEs is a result rather than a confound. The seed makes the set reproducible.
#
# Input   $GI_WORK/02_gene_sets/w22/biparental__gbM.bed      (02_gene_sets/03_build_sets.py)
#         $GI_DATA/references/W22v2.fa.fai                    (02_bisulfite_reference.sh)
# Output  $GI_WORK/04_emseq_methylation/metagene/shuffled_control.bed   (6,181 regions)
set -euo pipefail
source "$(dirname "$0")/../../config.sh"
gi_env gi_reads

GS=$GI_WORK/02_gene_sets/w22
WORK=$GI_WORK/04_emseq_methylation/metagene
mkdir -p "$WORK"

GFILE=$WORK/chr1_10.genome
# made by 02_bisulfite_reference.sh; made here when the mapping steps are skipped (README)
[ -s "$GI_DATA/references/W22v2.fa.fai" ] || samtools faidx "$GI_DATA/references/W22v2.fa"
awk 'BEGIN{OFS="\t"} $1 ~ /^chr([1-9]|10)$/ {print $1,$2}' "$GI_DATA/references/W22v2.fa.fai" > "$GFILE"

SHUF=$WORK/shuffled_control.bed
bedtools shuffle -i "$GS/biparental__gbM.bed" -g "$GFILE" -noOverlapping -seed 42 \
    | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"shuf_"NR,0,$6}' \
    | sort -k1,1 -k2,2n > "$SHUF"
echo "shuffled control: $(wc -l < "$SHUF") regions"
