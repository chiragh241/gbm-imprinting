#!/usr/bin/env bash
# Allele-resolved endosperm H3K27me3: align, assign reads to B73 or Mo17, count per gene body.
#
#   SRR748861 = B73 x Mo17 (B73 mother)   maternal = B73  = SNPsplit genome1
#   SRR748862 = Mo17 x B73 (Mo17 mother)  maternal = Mo17 = SNPsplit genome2
#
#   hisat2 --no-spliced-alignment --no-unal to the N-masked B73 v5 (HISAT2 2.2.1)
#   samtools view -q 30, coordinate sort
#   restricted to reads overlapping the informative windows, then SNPsplit (0.6.0)
#   featureCounts -s 0 -F SAF over gene bodies (Subread 2.0.8)
#
# A 36 bp read can be assigned to an allele only if it covers a variant, i.e. if it starts in
# [variant - 35, variant]; these windows, merged, are the "informative windows". Reads touching
# no window can never be assigned, so the BAM is restricted to them before SNPsplit, which is
# exact and makes SNPsplit finish in hours. Counting is over whole gene bodies, not exons:
# H3K27me3 is a broad domain and intronic reads carry the same allelic information.
#
# Input   $GI_DATA/reads/chip_zhang2014/<run>.fastq.gz                   (01_download_reads.sh)
#         $GI_WORK/08_rnaseq_allele/{index,snps,annot}                    (08_rnaseq_allele 02, 03)
# Output  $GI_WORK/10_chip_endosperm_allele/
#           annot/gene_bodies.saf, annot/informative_windows.bed
#           split/<run>.inf.genome{1,2}.bam
#           counts/<run>.{B73,Mo17}.txt, counts/<run>.total.txt
set -euo pipefail
source "$(dirname "$0")/../../config.sh"
gi_env gi_reads

M=$GI_WORK/08_rnaseq_allele
IDX=$M/index/B73v5_Nmasked
SNPF=$M/snps/snps.snpsplit.txt
SNPBED=$M/snps/snps.bed
W=$GI_WORK/10_chip_endosperm_allele
FQ=$GI_DATA/reads/chip_zhang2014
BAM=$W/bam; SPLIT=$W/split; CNT=$W/counts; ANN=$W/annot; LOG=$W/logs
mkdir -p "$BAM" "$SPLIT" "$CNT" "$ANN" "$LOG"

# gene-body SAF, chromosome names as in the genome ("1", not "chr1")
SAF=$ANN/gene_bodies.saf
if [ ! -s "$SAF" ]; then
    printf 'GeneID\tChr\tStart\tEnd\tStrand\n' > "$SAF"
    awk 'BEGIN{FS=OFS="\t"}{sub(/^chr/,"",$1); print $4,$1,$2+1,$3,$6}' "$M/annot/genes.bed" >> "$SAF"
fi

INF=$ANN/informative_windows.bed
if [ ! -s "$INF" ]; then
    awk 'BEGIN{OFS="\t"}{s=$2-35; if(s<0)s=0; print $1,s,$3}' "$SNPBED" \
        | sort -k1,1 -k2,2n | bedtools merge -i - > "$INF"
fi
echo "informative windows: $(wc -l < "$INF")  bp: $(awk '{s+=$3-$2}END{print s}' "$INF")"

for run in SRR748861 SRR748862; do
    echo ">>> $run $(date)"
    raw=$BAM/${run}.bam
    # BAMs are written under a temporary name and the allele split is marked when done, so an
    # interrupted run never leaves a truncated file for the next run to skip
    if [ ! -s "$raw" ]; then
        hisat2 -p "$GI_THREADS" --no-spliced-alignment --no-unal -U "$FQ/${run}.fastq.gz" -x "$IDX" \
            2> "$LOG/hisat2_${run}.log" \
            | samtools view -q 30 -bS -@ 2 - \
            | samtools sort -@ 4 -m 1G -O bam -o "$raw.part" -
        mv -f "$raw.part" "$raw"
        samtools index "$raw"
    fi
    [ -s "$CNT/${run}.total.txt" ] || featureCounts -s 0 -T 8 -F SAF -a "$SAF" \
        -o "$CNT/${run}.total.txt" "$raw" > "$LOG/fc_${run}_total.log" 2>&1
    inf=$BAM/${run}.inf.bam
    [ -s "$inf" ] || { samtools view -b -L "$INF" "$raw" > "$inf.part"; mv -f "$inf.part" "$inf"; }
    [ -f "$SPLIT/${run}.split_done" ] || { SNPsplit --snp_file "$SNPF" -o "$SPLIT" "$inf" \
        > "$LOG/snpsplit_${run}.log" 2>&1; touch "$SPLIT/${run}.split_done"; }
    for gn in 1 2; do
        lab=B73; [ "$gn" = 2 ] && lab=Mo17
        out=$CNT/${run}.${lab}.txt
        [ -s "$out" ] || featureCounts -s 0 -T 8 -F SAF -a "$SAF" -o "$out" \
            "$SPLIT/${run}.inf.genome${gn}.bam" > "$LOG/fc_${run}_${lab}.log" 2>&1
    done
done
ls -la "$CNT"
