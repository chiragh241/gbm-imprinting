#!/usr/bin/env bash
# One FASTQ pair per biological library: a replicate sequenced twice has its two runs
# concatenated. Merged files are recompressed as a single gzip member (zcat | pigz, not cat),
# so no reader can stop at the first member and silently drop the second run's reads.
#
# Input   $GI_DATA/reads/emseq_allele/<srr>_{1,2}.fastq.gz, data/allele/runs.tsv
# Output  $GI_WORK/05_allele_methylation/fastq/<library>_{1,2}.fastq.gz  (merged libraries)
#         $GI_WORK/05_allele_methylation/libs.tsv   library, tissue, direction, read 1, read 2
set -euo pipefail
source "$(dirname "$0")/../../config.sh"
gi_env gi_reads

FQ=$GI_DATA/reads/emseq_allele
LIST=$GI_SHIPPED/allele/runs.tsv
W=$GI_WORK/05_allele_methylation
MERGE=$W/fastq
LIBS=$W/libs.tsv
mkdir -p "$MERGE"

: > "$LIBS.tmp"
for lib in $(tail -n +2 "$LIST" | awk -F'\t' '!seen[$2]++ {print $2}'); do
    srrs=$(awk -F'\t' -v L="$lib" '$2==L{print $1}' "$LIST")
    tissue=$(awk -F'\t' -v L="$lib" '$2==L{print $3; exit}' "$LIST")
    direction=$(awk -F'\t' -v L="$lib" '$2==L{print $4; exit}' "$LIST")
    if [ "$(echo "$srrs" | wc -w)" -eq 1 ]; then
        r1=$FQ/${srrs}_1.fastq.gz; r2=$FQ/${srrs}_2.fastq.gz
    else
        r1=$MERGE/${lib}_1.fastq.gz; r2=$MERGE/${lib}_2.fastq.gz
        if [ ! -s "$r1" ] || [ ! -s "$r2" ]; then
            f1=""; f2=""
            for s in $srrs; do f1="$f1 $FQ/${s}_1.fastq.gz"; f2="$f2 $FQ/${s}_2.fastq.gz"; done
            echo "### $lib: merging $srrs"
            zcat $f1 | pigz -p "$GI_THREADS" > "$r1.part" && mv "$r1.part" "$r1"
            zcat $f2 | pigz -p "$GI_THREADS" > "$r2.part" && mv "$r2.part" "$r2"
        fi
    fi
    printf '%s\t%s\t%s\t%s\t%s\n' "$lib" "$tissue" "$direction" "$r1" "$r2" >> "$LIBS.tmp"
done
mv "$LIBS.tmp" "$LIBS"
cat "$LIBS"
