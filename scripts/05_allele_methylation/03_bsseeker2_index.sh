#!/usr/bin/env bash
# BS-Seeker2 bisulfite index of W22 v2 (bs_seeker2-build.py --aligner=bowtie2).
#
# The reference is not N-masked at B73/W22 variants: masking a C/T variant removes positions
# that bisulfite reads must be read through. Alleles are resolved from variant calls instead
# (08_snv.sh, 10_asm.sh).
#
# Input   $GI_DATA/references/W22v2.fa
# Output  $GI_WORK/05_allele_methylation/bsseeker_index/
set -euo pipefail
source "$(dirname "$0")/../../config.sh"
gi_env gi_bsseeker2

IDX=$GI_WORK/05_allele_methylation/bsseeker_index
mkdir -p "$IDX"
[ -f "$IDX/.build_done" ] && { echo "index already built"; exit 0; }
bs_seeker2-build.py -f "$GI_DATA/references/W22v2.fa" --aligner=bowtie2 -d "$IDX"
touch "$IDX/.build_done"
