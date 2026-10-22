#!/usr/bin/env bash
# EM-seq reads of W22 inbred tissues and demethylase mutants, BioProject PRJNA759188.
#
#   wild-type, mdr1 and dng102 (dng102-Q235) endosperm; wild-type 15-DAP embryo, developing
#   second leaf and premeiotic tassel; three biological replicates each, 24 runs.
#
# Input   data/emseq/runs.tsv   run accession, sample tag, size in the SRA run table
# Output  $GI_DATA/reads/emseq_w22/<tag>_1.fastq.gz, <tag>_2.fastq.gz
#
# Runs come from the SRA open-data mirror, are checked with vdb-validate and converted with
# fasterq-dump. vdb-validate rather than a size check: the run table's byte counts differ from
# the mirror's object sizes, so an exact size test would reject valid downloads. wget -c resumes
# a partial download; a finished run is skipped.
set -u
source "$(dirname "$0")/../../config.sh"
gi_env gi_reads

LIST=$GI_SHIPPED/emseq/runs.tsv
FQ=$GI_DATA/reads/emseq_w22
SRA=$GI_WORK/04_emseq_methylation/sra
TMP=$GI_WORK/04_emseq_methylation/tmp
mkdir -p "$FQ" "$SRA" "$TMP"

while IFS=$'\t' read -r run tag bytes; do
    [ -z "${run:-}" ] && continue
    if [ -s "$FQ/${tag}_1.fastq.gz" ] && [ -s "$FQ/${tag}_2.fastq.gz" ]; then
        echo "### skip $run ($tag), present"; continue
    fi
    echo "### $run -> $tag  $(date)"
    wget -q -c "https://sra-pub-run-odp.s3.amazonaws.com/sra/$run/$run" -O "$SRA/$run" \
        || { echo "!!! download failed $run"; exit 1; }
    vdb-validate "$SRA/$run" || { echo "!!! vdb-validate failed $run"; exit 1; }
    fasterq-dump "$SRA/$run" --split-files --force --threads "$GI_THREADS" --temp "$TMP" \
        --outdir "$FQ" || { echo "!!! fasterq-dump failed $run"; exit 1; }
    [ -s "$FQ/${run}_1.fastq" ] && [ -s "$FQ/${run}_2.fastq" ] \
        || { echo "!!! $run is not paired"; exit 1; }
    # compressed under a temporary name, so an interrupted run never leaves a truncated file
    # under the final name for the next run to skip
    for m in 1 2; do
        pigz -p "$GI_THREADS" -c "$FQ/${run}_$m.fastq" > "$FQ/${tag}_$m.fastq.gz.part" \
            && mv -f "$FQ/${tag}_$m.fastq.gz.part" "$FQ/${tag}_$m.fastq.gz" \
            || { echo "!!! compression failed $run"; exit 1; }
    done
    rm -f "$FQ/${run}_1.fastq" "$FQ/${run}_2.fastq" "$SRA/$run"
    rm -rf "${TMP:?}"/*
done < <(tail -n +2 "$LIST")
echo "### fastq pairs present: $(ls "$FQ"/*_1.fastq.gz 2>/dev/null | wc -l) of 24"
