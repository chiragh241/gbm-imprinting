#!/usr/bin/env bash
# Raw runs (data/rice/runs.tsv) as .sra files from NCBI's copy on AWS, md5-checked against NCBI's
# locate service (sra_locate.py): the eight 5 DAP F1 endosperm WGBS libraries of Tonosaki et al.
# 2024 (DRA015356; Nipponbare x Kasalath, Kasalath x Nipponbare, Nipponbare x 93-11, 93-11 x
# Nipponbare, two replicates each) and the Kasalath and 93-11 resequencing runs they list
# (DRR013715/6, SRR933668/9). About 210 GB. The .sra files keep the base qualities the trimming
# needs and are read directly by fastq-dump, so no FASTQ is written.
#
# usage: 02_download_reads.sh [reseq|wgbs|all] [parallel downloads, default 3]
#        02_download_reads.sh RUN [RUN ...]
# Output $GI_DATA/rice/sra/<run>.sra, $GI_DATA/rice/sra/sra_md5.tsv (run, md5, size)
#
# An interrupted download keeps its .part and resumes from it; a .part is deleted only when it is
# complete in size and fails the md5.
set -euo pipefail
source "$(dirname "$0")/rice_env.sh"
gi_env gi_rice
mkdir -p "$SRA"
HERE=$(cd "$(dirname "$0")" && pwd)

fetch() {
    local run=$1 out=$SRA/$1.sra info url md5 size have
    [ -s "$out" ] && { echo "$run: present"; return 0; }
    info=$(python "$HERE/sra_locate.py" "$run")
    url=$(cut -f1 <<< "$info"); md5=$(cut -f2 <<< "$info"); size=$(cut -f3 <<< "$info")
    echo "$(date +%T) $run: $size bytes from $url"
    if ! curl -fsSL --retry 20 --retry-delay 30 -C - -o "$out.part" "$url"; then
        echo "!!! $run: download interrupted at $(stat -c%s "$out.part" 2>/dev/null || echo 0) bytes, .part kept"
        return 1
    fi
    have=$(stat -c%s "$out.part")
    if [ "$have" != "$size" ]; then
        echo "!!! $run: $have of $size bytes, .part kept"; return 1
    fi
    if ! echo "$md5  $out.part" | md5sum -c --quiet -; then
        echo "!!! $run: md5 mismatch, .part removed"; rm -f "$out.part"; return 1
    fi
    mv "$out.part" "$out"
    printf "%s\t%s\t%s\n" "$run" "$md5" "$size" >> "$SRA/sra_md5.tsv"
    echo "$(date +%T) $run: done"
}
export -f fetch
export SRA HERE

case "${1:-all}" in
    reseq|wgbs|all)
        awk -F'\t' -v t="${1:-all}" 'NR > 1 && (t == "all" || $2 == t) {print $1}' "$RUNS" \
            | xargs -P "${2:-3}" -I{} bash -c 'fetch {}' ;;
    *)
        printf "%s\n" "$@" | xargs -P 3 -I{} bash -c 'fetch {}' ;;
esac
