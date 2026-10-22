#!/usr/bin/env bash
# B73/W22 single-nucleotide variants from the pooled endosperm ATCGmaps (cgmaptools snv).
#
# Sites are kept at read depth 15 to 70: shallower genotypes are noise, deeper sites are mostly
# collapsed repeats that would mis-phase reads. C/T and G/A variants are kept on purpose: at a
# C/T site a read still showing C must be the unconverted C allele, so they stay informative.
#
# Input   $GI_WORK/05_allele_methylation/atcgmap/<sequence>.ATCGmap.gz   (07_call_methylation.sh)
# Output  $GI_WORK/05_allele_methylation/snv/endosperm_snv.vcf   (5,158,268 variants)
#         per-sequence VCFs are concatenated in file-name order (chr1, chr10, chr2, ... unmapped)
set -u
source "$(dirname "$0")/../../config.sh"
gi_env gi_cgmaptools

W=$GI_WORK/05_allele_methylation
ATCG=$W/atcgmap
OUT=$W/snv
VCF=$OUT/endosperm_snv.vcf
MINCOV=15; MAXCOV=70
MAXJ=${GI_SNV_JOBS:-4}
mkdir -p "$OUT"

mapfile -t MAPS < <(ls "$ATCG"/*.ATCGmap.gz)
for m in "${MAPS[@]}"; do
    chr=$(basename "$m" .ATCGmap.gz)
    [ -s "$OUT/${chr}.vcf" ] && continue
    while [ "$(jobs -rp | wc -l)" -ge "$MAXJ" ]; do wait -n; done
    (
        raw=$OUT/${chr}.raw.vcf
        # the tabular result on stdout is not needed, only the VCF written with -v
        zcat "$m" | cgmaptools snv -v "$raw" > /dev/null 2> "$OUT/${chr}.snv.log"
        prc=("${PIPESTATUS[@]}")
        if [ "${prc[0]}" -eq 0 ] && [ "${prc[1]}" -eq 0 ] && [ -s "$raw" ]; then
            # depth is in INFO as "NS=1:DP=86"
            # under a temporary name, so an interrupted run leaves no truncated VCF to skip
            awk -v lo="$MINCOV" -v hi="$MAXCOV" '
                /^#/ { print; next }
                { dp = $8; sub(/.*DP=/, "", dp); sub(/[^0-9].*/, "", dp)
                  if (dp+0 >= lo && dp+0 <= hi) print }' "$raw" > "$OUT/${chr}.vcf.part" \
                && mv -f "$OUT/${chr}.vcf.part" "$OUT/${chr}.vcf"
            rm -f "$raw"
        else
            echo "!!! $chr failed (see $OUT/${chr}.snv.log)"
        fi
    ) &
done
wait

for m in "${MAPS[@]}"; do
    [ -s "$OUT/$(basename "$m" .ATCGmap.gz).vcf" ] || { echo "!!! missing $m"; exit 1; }
done
first=1; : > "$VCF"
for m in "${MAPS[@]}"; do
    f=$OUT/$(basename "$m" .ATCGmap.gz).vcf
    if [ "$first" = 1 ]; then grep '^#' "$f" >> "$VCF"; first=0; fi
    grep -v '^#' "$f" >> "$VCF"
done
echo "$(grep -vc '^#' "$VCF") variants in $VCF"
