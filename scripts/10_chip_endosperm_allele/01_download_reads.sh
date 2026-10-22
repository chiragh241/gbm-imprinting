#!/usr/bin/env bash
# Allele-specific endosperm H3K27me3 ChIP-seq of Zhang et al. 2014 (SRP011991), 36 bp single-end:
#   SRR748861  B73 x Mo17 endosperm, 12 DAP (B73 mother)   12.8 GB
#   SRR748862  Mo17 x B73 endosperm, 12 DAP (Mo17 mother)  15.4 GB
# No input, mock or IgG library exists for these samples.
#
# Output  $GI_DATA/reads/chip_zhang2014/<run>.fastq.gz
#
# ENA's HTTPS server tends to drop long transfers, so each file is resumed alternately over
# HTTPS and FTP until its size matches ENA's exactly, and is then checked against ENA's md5.
set -uo pipefail
source "$(dirname "$0")/../../config.sh"
FQ=$GI_DATA/reads/chip_zhang2014
mkdir -p "$FQ"

declare -A B=( [SRR748861]=12815732651 [SRR748862]=15373211886 )
declare -A M=( [SRR748861]=e2578198ae0930bdbb8ee49bd4f97e2f [SRR748862]=ec675d8d59cd086ce06886337064d80d )
for run in SRR748861 SRR748862; do
    out=$FQ/${run}.fastq.gz
    want=${B[$run]}
    path="SRR748/$run/${run}.fastq.gz"
    for ((i = 1; i <= 300; i++)); do
        got=$(stat -c%s "$out" 2>/dev/null || echo 0)
        [ "$got" -eq "$want" ] && break
        [ "$got" -gt "$want" ] && { rm -f "$out"; got=0; }
        if (( i % 2 == 1 )); then url="https://ftp.sra.ebi.ac.uk/vol1/fastq/$path"
        else                      url="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$path"; fi
        wget -c -q --tries=0 --timeout=120 --read-timeout=120 --waitretry=15 \
             --retry-connrefused -O "$out" "$url"
        sleep 5
    done
    [ "$(stat -c%s "$out")" -eq "$want" ] || { echo "!!! $run incomplete"; exit 1; }
    [ "$(md5sum < "$out" | cut -c1-32)" = "${M[$run]}" ] \
        || { echo "!!! $run does not match ENA's md5; delete $out and run again"; exit 1; }
    echo "$run complete, md5 ok"
done
