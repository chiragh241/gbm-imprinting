#!/usr/bin/env bash
# Allele-specific methylation of each endosperm library at the kept variants (cgmaptools asm).
#
#   cgmaptools asm -m ass -d 2 -t <CG|CHG|CHH>
#     -m ass  allele-specific site calling; -d 2  at least two reads per allele at a site
#
# Each library is phased from its own deduplicated BAM. Maternal and paternal are assigned
# downstream from the cross direction (11_gene_summary.py), not from the genotype.
#
# Input   $GI_WORK/05_allele_methylation/dedup/<library>.dedup.bam     (05_sort_dedup.sh)
#         $GI_WORK/05_allele_methylation/snv/endosperm_snv.asm.vcf      (09_asm_snv_subset.py)
#         $GI_DATA/references/W22v2.fa (+ .fai)
# Output  $GI_WORK/05_allele_methylation/asm/<library>.<context>.asm.tsv   (18 tables)
#         These tables are deposited on Zenodo (see README).
set -u
source "$(dirname "$0")/../../config.sh"
gi_env gi_cgmaptools

W=$GI_WORK/05_allele_methylation
VCF=$W/snv/endosperm_snv.asm.vcf
REF=$GI_DATA/references/W22v2.fa
OUT=$W/asm
MAXJ=${GI_ASM_JOBS:-6}
mkdir -p "$OUT"
command -v R >/dev/null && command -v samtools >/dev/null \
    || { echo "!!! cgmaptools asm needs R and samtools on PATH"; exit 1; }

for l in $(awk -F'\t' '$2=="endosperm"{print $1}' "$W/libs.tsv"); do
    for ctx in CG CHG CHH; do
        [ -f "$OUT/${l}.${ctx}.asm_done" ] && continue
        while [ "$(jobs -rp | wc -l)" -ge "$MAXJ" ]; do wait -n; done
        (
            o=$OUT/${l}.${ctx}.asm.tsv
            if cgmaptools asm -r "$REF" -b "$W/dedup/${l}.dedup.bam" -l "$VCF" \
                   -t "$ctx" -m ass -d 2 > "$o" 2> "$OUT/${l}.${ctx}.log" && [ -s "$o" ]; then
                touch "$OUT/${l}.${ctx}.asm_done"
            else
                echo "!!! $l $ctx failed (see $OUT/${l}.${ctx}.log)"
            fi
        ) &
    done
done
wait
n=$(ls "$OUT"/*.asm_done 2>/dev/null | wc -l)
echo "$n of 18 library x context tables"
[ "$n" -eq 18 ]
