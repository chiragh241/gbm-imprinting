#!/usr/bin/env bash
# Tonosaki et al.'s 5 DAP endosperm WGBS libraries, each from .sra to per-allele CG methylation.
#
#   1 read pairs streamed from the .sra (fastq-dump) through cutadapt with Trim Galore's defaults
#     ($TRIM_ARGS; Tonosaki used Cutadapt 1.17 and Trimmomatic 0.36) into two gzipped files
#   2 Bismark (Bowtie2) against IRGSP-1.0 with this cross's SNPs set to N (04), $BISMARK_OPTS,
#     temporary files gzipped. --multicore stays at 5: Bismark merges its per-process outputs in
#     process order, and deduplication keeps the first copy of a duplicate it meets, so another
#     process count would keep other copies of some duplicated pairs
#   3 deduplicate_bismark, one chromosome at a time (see step 3)
#   4 SNPsplit --bisulfite: read pairs carrying Nipponbare bases -> genome1, the other parent's
#     bases -> genome2 (a pair is assigned when at least one read covers a SNP and none conflicts;
#     Tonosaki: "Both paired reads were also kept if one pair of reads covered an SNP")
#   5 bismark_methylation_extractor --no_overlap on each allele file -> per-CG coverage file
#   6 bisulfite conversion from the chloroplast reads, which are unmethylated
#
# Each step leaves a mark, so a rerun resumes after the last finished step. Trimmed reads and the
# extractor's per-read call files are deleted once used; the BAMs and coverage files are kept.
# Here Bismark took 15 to 17 hours per library and steps 3 to 6 another 7 to 9, about 8 days for
# the eight, with about 120 GB of BAMs per library; the eight libraries' coverage files and
# reports are deposited on Zenodo (see README).
#
# usage: 05_map_libraries.sh [RUN ...]   (default: every WGBS run in data/rice/runs.tsv)
# Input   $GI_DATA/rice/sra/<run>.sra (02), $GI_WORK/13_rice/nmask_*/ (04),
#         $GI_WORK/13_rice/snps/snps_*.txt (03)
# Output  $GI_WORK/13_rice/bam/<run>*, $GI_WORK/13_rice/meth/<run>.genome{1,2}.bismark.cov.gz,
#         reports in $GI_WORK/13_rice/qc/
set -euo pipefail
source "$(dirname "$0")/rice_env.sh"
gi_env gi_rice
T=$RW/trim; B=$RW/bam; M=$RW/meth; Q=$RW/qc
mkdir -p "$T" "$B" "$M" "$Q"
SPLIT_AWK=$(command -v mawk || command -v awk)      # mawk splits about twice as fast

map_library() {
    local run=$1 alt G R1 R2 BAM DEDUP STEM C g
    alt=$(awk -F'\t' -v r="$run" '$1 == r && $2 == "wgbs" {print $3}' "$RUNS")
    [ -n "$alt" ] || { echo "$run is not a WGBS run in $RUNS"; return 1; }
    G=$RW/nmask_$alt
    mkdir -p "$TMPDIR/$run"
    mark() { touch "$B/$run.$1.ok"; echo "$(date '+%F %T') $run: $1 done"; }
    done_() { [ -e "$B/$run.$1.ok" ]; }

    # 1 reads, trimmed
    R1=$T/${run}_1_val_1.fq.gz
    R2=$T/${run}_2_val_2.fq.gz
    if ! done_ trim; then
        local qb
        qb=$(qual_base "$run")
        echo "$(date '+%F %T') $run: Phred+$qb"
        # shellcheck disable=SC2086
        fastq-dump --split-spot --stdout "${SPOTS[@]}" "$SRA/$run.sra" 2> "$T/$run.fastq-dump.log" \
            | cutadapt -j "$GI_THREADS" --quality-base "$qb" $TRIM_ARGS --interleaved \
                  -o "$T/${run}_1.part.fq.gz" -p "$T/${run}_2.part.fq.gz" - > "$T/${run}_cutadapt_report.txt"
        mv "$T/${run}_1.part.fq.gz" "$R1"
        mv "$T/${run}_2.part.fq.gz" "$R2"
        mark trim
    fi
    BAM=$B/${run}_1_val_1_bismark_bt2_pe.bam
    DEDUP=$B/${run}_1_val_1_bismark_bt2_pe.deduplicated.bam
    STEM=$B/${run}_1_val_1_bismark_bt2_pe.deduplicated

    # 2 Bismark
    if ! done_ bismark; then
        rm -f "$BAM"
        # shellcheck disable=SC2086
        bismark --genome "$G" $BISMARK_OPTS --multicore 5 --temp_dir "$TMPDIR/$run" --gzip \
            -o "$B" -1 "$R1" -2 "$R2"
        samtools quickcheck "$BAM"
        mark bismark
        rm -f "$R1" "$R2"
    fi

    # 3 duplicates, one chromosome at a time. deduplicate_bismark keeps every unique fragment of
    #   its input in memory, which on a whole 190M-pair library grows past 12 GB within minutes.
    #   Both reads of a pair, and all its duplicates, share a chromosome, and splitting keeps the
    #   read order, so the pairs removed are exactly those of one whole-file run. Two chromosomes
    #   at a time (~3 GB each).
    if ! done_ dedup; then
        C=$B/${run}_bychr
        rm -rf "$C"; mkdir -p "$C"
        samtools view -h -@ 4 "$BAM" | "$SPLIT_AWK" -v d="$C" '
            /^@/ { h = h $0 "\n"; next }
            { if (!($3 in out)) { out[$3] = "samtools view -b -o " d "/" $3 ".bam -"; printf "%s", h | out[$3] }
              print | out[$3] }
            END { for (c in out) close(out[c]) }'
        printf "%s\n" "$C"/*.bam | xargs -P 2 -I{} deduplicate_bismark --paired --bam --output_dir "$C" {} \
            > "$C/deduplicate_bismark.log" 2>&1
        # shellcheck disable=SC2046
        samtools cat -o "$DEDUP.part" $(find "$C" -name '*.deduplicated.bam' | sort)
        mv "$DEDUP.part" "$DEDUP"
        samtools quickcheck "$DEDUP"
        { cat "$C"/*.deduplication_report.txt
          awk -F'\t' '/^Total number of alignments analysed/ {n += $2}
                      /^Total number duplicated alignments removed/ {split($2, a, " "); d += a[1]}
                      END {printf "All chromosomes: %d pairs analysed, %d duplicates removed (%.2f%%), %d left\n",
                           n, d, 100 * d / n, n - d}' "$C"/*.deduplication_report.txt
        } > "$B/${run}_1_val_1_bismark_bt2_pe.deduplication_report.txt"
        tail -n 1 "$B/${run}_1_val_1_bismark_bt2_pe.deduplication_report.txt"
        rm -rf "$C"
        mark dedup
    fi

    # 4 allele split
    if ! done_ snpsplit; then
        SNPsplit --paired --bisulfite --snp_file "$RW/snps/snps_$alt.txt" -o "$B" "$DEDUP"
        samtools quickcheck "$STEM.genome1.bam" "$STEM.genome2.bam"
        mark snpsplit
    fi

    # 5 per-allele CG calls
    for g in genome1 genome2; do
        if ! done_ "extract_$g"; then
            bismark_methylation_extractor --paired-end --no_overlap --comprehensive --gzip --bedGraph \
                --buffer_size 8G --parallel 4 -o "$M" "$STEM.$g.bam"
            mv "$M/${run}_1_val_1_bismark_bt2_pe.deduplicated.$g.bismark.cov.gz" "$M/$run.$g.bismark.cov.gz"
            rm -f "$M"/C*_context_"${run}"_1_val_1_bismark_bt2_pe.deduplicated."$g".txt.gz
            mark "extract_$g"
        fi
    done

    # 6 bisulfite conversion from the chloroplast (unmethylated), as Tonosaki estimated it
    if ! done_ chrC; then
        samtools view -h "$DEDUP" | awk '$1 ~ /^@/ || $3 == "ChrC"' | samtools view -b -o "$B/$run.chrC.bam" -
        bismark_methylation_extractor --paired-end --no_overlap --comprehensive --gzip -o "$M/chrC" \
            "$B/$run.chrC.bam"
        rm -f "$M/chrC"/C*_context_"$run".chrC.txt.gz
        mark chrC
    fi

    # reports, for 09_qc_table.py
    cp "$B/${run}_1_val_1_bismark_bt2_PE_report.txt" "$B/${run}"*deduplication_report.txt \
       "$B/${run}"*SNPsplit_report.txt "$B/${run}"*SNPsplit_sort.txt "$M/${run}"*M-bias.txt \
       "$M/${run}"*splitting_report.txt "$Q/" 2>/dev/null || true
    cp "$M/chrC/$run.chrC_splitting_report.txt" "$T/${run}_cutadapt_report.txt" "$Q/" 2>/dev/null || true
    touch "$B/$run.done"
    echo "$(date '+%F %T') $run: finished"
}

runs=("$@")
[ ${#runs[@]} -gt 0 ] || mapfile -t runs < <(awk -F'\t' '$2 == "wgbs" {print $1}' "$RUNS")
for run in "${runs[@]}"; do
    [ -e "$B/$run.done" ] && { echo "### skip $run"; continue; }
    echo "=== $(date '+%F %T') $run"
    map_library "$run"
done
