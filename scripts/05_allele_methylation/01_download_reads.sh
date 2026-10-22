#!/usr/bin/env bash
# EM-seq reads of reciprocal B73 x W22 and W22 x B73 F1 mature endosperm, BioProject PRJNA759188,
# three biological replicates per direction; four replicates were sequenced twice (orig + reseq).
#
# Input   data/allele/runs.tsv   srr, library, tissue, direction (mother x father), rep, role
# Output  $GI_DATA/reads/emseq_allele/<srr>_1.fastq.gz, <srr>_2.fastq.gz
set -u
source "$(dirname "$0")/../../config.sh"
gi_env gi_reads

FQ=$GI_DATA/reads/emseq_allele
SRA=$GI_WORK/05_allele_methylation/sra
TMP=$GI_WORK/05_allele_methylation/tmp
mkdir -p "$FQ" "$SRA" "$TMP"

while IFS=$'\t' read -r run library tissue direction rep role; do
    [ -z "${run:-}" ] && continue
    if [ -s "$FQ/${run}_1.fastq.gz" ] && [ -s "$FQ/${run}_2.fastq.gz" ]; then
        echo "### skip $run, present"; continue
    fi
    echo "### $run ($library $role)  $(date)"
    wget -q -c "https://sra-pub-run-odp.s3.amazonaws.com/sra/$run/$run" -O "$SRA/$run" \
        || { echo "!!! download failed $run"; exit 1; }
    vdb-validate "$SRA/$run" || { echo "!!! vdb-validate failed $run"; exit 1; }
    fasterq-dump "$SRA/$run" --split-files --force --threads "$GI_THREADS" --temp "$TMP" \
        --outdir "$FQ" || { echo "!!! fasterq-dump failed $run"; exit 1; }
    [ -s "$FQ/${run}_1.fastq" ] && [ -s "$FQ/${run}_2.fastq" ] || { echo "!!! $run is not paired"; exit 1; }
    # compressed under a temporary name, so an interrupted run never leaves a truncated file
    # under the final name for the next run to skip
    for m in 1 2; do
        pigz -p "$GI_THREADS" -c "$FQ/${run}_$m.fastq" > "$FQ/${run}_$m.fastq.gz.part" \
            && mv -f "$FQ/${run}_$m.fastq.gz.part" "$FQ/${run}_$m.fastq.gz" \
            || { echo "!!! compression failed $run"; exit 1; }
    done
    rm -f "$FQ/${run}_1.fastq" "$FQ/${run}_2.fastq" "$SRA/$run"; rm -rf "${TMP:?}"/*
done < <(tail -n +2 "$GI_SHIPPED/allele/runs.tsv")
echo "### fastq pairs present: $(ls "$FQ"/*_1.fastq.gz 2>/dev/null | wc -l) of 10"
