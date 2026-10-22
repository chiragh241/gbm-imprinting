#!/usr/bin/env bash
# One Bismark genome per cross: IRGSP-1.0 with that cross's SNP positions set to N (nmask.py), as
# Tonosaki et al. did ("reference Nipponbare genome (IRGSP-1.0) that N-masked on SNP position")
# and as SNPsplit expects, so reads from either parent align without favouring the Nipponbare base.
#
# Input   $GI_DATA/rice/ref/genome.fa (01), $GI_WORK/13_rice/snps/snps_{KAS,9311}.txt (03)
# Output  $GI_WORK/13_rice/nmask_{KAS,9311}/genome_nmask_*.fa and its Bisulfite_Genome/
#         $GI_WORK/13_rice/snps/nmask_{KAS,9311}.tsv   SNPs masked per sequence
set -euo pipefail
source "$(dirname "$0")/rice_env.sh"
gi_env gi_rice
HERE=$(cd "$(dirname "$0")" && pwd)
for acc in KAS 9311; do
    G=$RW/nmask_$acc
    mkdir -p "$G"
    if [ ! -s "$G/genome_nmask_$acc.fa" ]; then
        python "$HERE/nmask.py" "$REF/genome.fa" "$RW/snps/snps_$acc.txt" "$G/genome_nmask_$acc.fa.part" \
            > "$RW/snps/nmask_$acc.tsv"
        mv "$G/genome_nmask_$acc.fa.part" "$G/genome_nmask_$acc.fa"
    fi
    if [ ! -e "$G/Bisulfite_Genome/.done" ]; then
        rm -rf "$G/Bisulfite_Genome"
        bismark_genome_preparation --bowtie2 --parallel 6 "$G" > "$G/genome_preparation.log" 2>&1
        touch "$G/Bisulfite_Genome/.done"
    fi
done
