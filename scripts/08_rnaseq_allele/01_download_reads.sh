#!/usr/bin/env bash
# Endosperm RNA-seq of reciprocal B73 x Mo17 crosses, the B73 parent wild type or the
# zmros1a zmros1b (dng103 mdr1) double mutant (Xu et al. 2022, PRJNA739488): 12 runs, three
# replicates of each of four crosses (data/rnaseq/samplesheet.tsv). Paired-end 2 x 150 bp.
#
# Input   data/rnaseq/fastq_md5.tsv   ENA FASTQ URL and md5 for each of the 24 files
# Output  $GI_DATA/reads/rnaseq_xu2022/<run>_{1,2}.fastq.gz   (about 53 GB)
#
# Every file is checked against ENA's md5; a file that already matches is not fetched again.
set -u
source "$(dirname "$0")/../../config.sh"
FQ=$GI_DATA/reads/rnaseq_xu2022
mkdir -p "$FQ"

bad=0
while IFS=$'\t' read -r url want; do
    [ -n "${url:-}" ] || continue
    f=$FQ/$(basename "$url")
    if [ -s "$f" ] && [ "$(md5sum "$f" | awk '{print $1}')" = "$want" ]; then
        echo "  ok   $(basename "$f")"; continue
    fi
    wget -q -c -P "$FQ" "$url" || { echo "  FAILED $url"; bad=1; continue; }
    if [ "$(md5sum "$f" | awk '{print $1}')" = "$want" ]; then
        echo "  got  $(basename "$f")"
    else
        echo "  BAD md5 $(basename "$f")"; bad=1
    fi
done < "$GI_SHIPPED/rnaseq/fastq_md5.tsv"
exit $bad
