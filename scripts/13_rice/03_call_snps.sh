#!/usr/bin/env bash
# Kasalath and 93-11 SNPs against IRGSP-1.0, from the resequencing runs Tonosaki et al. list in
# Supplementary Data 1 (Kasalath DRR013715/6, 93-11 SRR933668/9).
#
# Their pipeline: Cutadapt and Trimmomatic, Bowtie2 to IRGSP-1.0, Picard duplicate removal, GATK
# 3.8 (realignment, recalibration, calling); positions kept where at least 90% of five or more
# reads carry the SNP; exonic SNPs then corrected with RNA-seq of selfed Kasalath and 93-11
# endosperm. Here: cutadapt with Trim Galore's defaults, Bowtie2, samtools markdup, bcftools
# mpileup/call on both accessions jointly; kept per accession: biallelic SNPs, homozygous
# alternative, depth >= 5 and the alternative base on >= 90% of reads. The RNA-seq correction is
# not repeated (those selfed endosperm runs are not among their listed accessions). Their counts:
# NIP-KAS 1,664,378 and NIP-9311 1,753,321 SNPs; here 1,883,842 and 1,998,737.
#
# Input   $GI_DATA/rice/sra/<run>.sra (02), $GI_DATA/rice/ref/genome.fa (01)
# Output  $GI_WORK/13_rice/snps/snps_{KAS,9311}.txt   SNPsplit format (ID, Chr, Position, SNP
#                                                    value, Ref/SNP)
#         $GI_WORK/13_rice/snps/snp_counts.tsv
set -euo pipefail
source "$(dirname "$0")/rice_env.sh"
gi_env gi_rice
S=$RW/reseq
O=$RW/snps
mkdir -p "$S/calls" "$RW/bt2" "$O"

# 1 Bowtie2 index of the unmasked reference
if [ ! -s "$RW/bt2/genome.rev.2.bt2" ]; then
    bowtie2-build --threads "$GI_THREADS" "$REF/genome.fa" "$RW/bt2/genome" > "$S/bt2_build.log" 2>&1
fi

# 2 per run: read pairs streamed from the .sra through cutadapt into Bowtie2 (no FASTQ on disk),
#   then fixmate and sort; per accession: merge, drop duplicates
for acc in KAS 9311; do
    runs=$(awk -F'\t' -v a="$acc" '$2 == "reseq" && $3 == a {print $1}' "$RUNS")
    for run in $runs; do
        bam=$S/$run.bam
        [ -s "$bam" ] && continue
        qb=$(qual_base "$run")
        p64=""; [ "$qb" = 64 ] && p64="--phred64"
        echo "$(date '+%F %T') $run: Phred+$qb"
        # shellcheck disable=SC2086
        fastq-dump --split-spot --stdout "${SPOTS[@]}" "$SRA/$run.sra" 2> "$S/$run.fastq-dump.log" \
            | cutadapt -j 4 --quality-base "$qb" $TRIM_ARGS --interleaved -o - - 2> "$S/$run.cutadapt.log" \
            | bowtie2 -p "$GI_THREADS" $p64 --rg-id "$run" --rg "SM:$acc" --rg PL:ILLUMINA \
                  -x "$RW/bt2/genome" --interleaved - 2> "$S/$run.bowtie2.log" \
            | samtools fixmate -m -@ 2 - - \
            | samtools sort -@ 4 -m 1500M -T "$TMPDIR/$run" -o "$bam.part" -
        mv "$bam.part" "$bam"
    done
    md=$S/$acc.markdup.bam
    if [ ! -s "$md" ]; then
        # shellcheck disable=SC2046
        samtools merge -@ 4 -f -o "$S/$acc.merged.bam" $(for r in $runs; do echo "$S/$r.bam"; done)
        samtools markdup -@ 4 -r -s "$S/$acc.merged.bam" "$md.part" 2> "$S/$acc.markdup.log"
        mv "$md.part" "$md"
        samtools index "$md"
        rm -f "$S/$acc.merged.bam"
    fi
done

# 3 joint calling, one job per chromosome (organelles left out)
call() {
    local c=$1 out=$S/calls/$1.vcf.gz
    [ -s "$out" ] && return 0
    bcftools mpileup -f "$REF/genome.fa" -r "$c" -q 20 -Q 20 -d 1000 -a FORMAT/AD,FORMAT/DP -Ou \
            "$S/KAS.markdup.bam" "$S/9311.markdup.bam" \
        | bcftools call -m -v -Oz -o "$out.part"
    mv "$out.part" "$out"
}
export -f call
export S REF
chroms=$(cut -f1 "$REF/genome.fa.fai" | grep -v -x -e ChrC -e ChrM)
printf "%s\n" $chroms | xargs -P "$GI_THREADS" -I{} bash -c 'call {}'
if [ ! -s "$S/raw.vcf.gz" ]; then
    # shellcheck disable=SC2046
    bcftools concat -Oz -o "$S/raw.vcf.gz.part" $(for c in $chroms; do echo "$S/calls/$c.vcf.gz"; done)
    mv "$S/raw.vcf.gz.part" "$S/raw.vcf.gz"
    bcftools index -t "$S/raw.vcf.gz"
fi

# 4 per accession: homozygous alternative SNPs on >= 5 reads, >= 90% of them alternative
printf "accession\tsnps\n" > "$O/snp_counts.tsv"
for acc in KAS 9311; do
    out=$O/snps_$acc.txt
    if [ ! -s "$out" ]; then
        bcftools view -s "$acc" -a -Ou "$S/raw.vcf.gz" \
            | bcftools view -v snps -m2 -M2 \
                  -i 'GT="AA" && FMT/DP>=5 && FMT/AD[0:1] >= 0.9*(FMT/AD[0:0]+FMT/AD[0:1])' -Ou \
            | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
            | awk 'BEGIN {OFS = "\t"; print "ID", "Chr", "Position", "SNP value", "Ref/SNP"}
                   {print NR, $1, $2, 1, $3 "/" $4}' > "$out.part"
        mv "$out.part" "$out"
    fi
    printf "%s\t%s\n" "$acc" "$(($(wc -l < "$out") - 1))" >> "$O/snp_counts.tsv"
done
cat "$O/snp_counts.tsv"
