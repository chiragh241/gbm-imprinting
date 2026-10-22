#!/usr/bin/env bash
# Rice reference and the published tables of Tonosaki et al. 2024 (Nature Plants 10:1231).
#
#   IRGSP-1.0 pseudomolecules and MSU Release 7 gene models, as served by the Rice Genome
#   Annotation Project (rice.uga.edu). Tonosaki et al. mapped to IRGSP-1.0 and name genes by MSU
#   locus (LOC_Os...). The Nipponbare chloroplast (NC_001320.1) and mitochondrion (NC_011033.1)
#   from NCBI are added as ChrC and ChrM, so organelle reads map to their own genome rather than
#   to nuclear insertions, and ChrC gives the bisulfite conversion rate (as Tonosaki estimated it).
#   locus_brief_info flags the MSU7 loci annotated as transposable element related (is_TE).
#
#   Supplementary Data 1 (sequencing runs and per-library numbers) and 3 (imprint class of every
#   MEG and PEG at 3, 5 and 7 DAP), saved under the names the scripts read, md5-checked.
#
# Output  $GI_DATA/rice/ref/genome.fa (+ .fai), osa1_r7.all_models.gff3.gz,
#         osa1_r7.locus_brief_info.txt.gz, reference_sha256.txt
#         $GI_DATA/rice/supplementary/tonosaki2024_SupplementaryData{1_NGS,3_imprint_classes}.xlsx
set -euo pipefail
source "$(dirname "$0")/rice_env.sh"
gi_env gi_rice
mkdir -p "$REF" "$SUPP"

U=https://rice.uga.edu/osa1r7_download
for f in osa1_r7.asm.fa.gz osa1_r7.all_models.gff3.gz osa1_r7.locus_brief_info.txt.gz; do
    gi_fetch "$U/$f" "$REF/$f"
done
E="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text"
gi_fetch "$E&id=NC_001320.1" "$REF/NC_001320.1.fa"    # Oryza sativa Japonica Group chloroplast
gi_fetch "$E&id=NC_011033.1" "$REF/NC_011033.1.fa"    # Oryza sativa Japonica Group mitochondrion

if [ ! -s "$REF/genome.fa" ]; then
    { zcat "$REF/osa1_r7.asm.fa.gz" | sed 's/^>\([^ ]*\).*/>\1/'
      sed '1s/.*/>ChrC/' "$REF/NC_001320.1.fa"
      sed '1s/.*/>ChrM/' "$REF/NC_011033.1.fa"; } > "$REF/genome.fa.part"
    mv "$REF/genome.fa.part" "$REF/genome.fa"
fi
[ -s "$REF/genome.fa.fai" ] || samtools faidx "$REF/genome.fa"
( cd "$REF" && sha256sum osa1_r7.asm.fa.gz osa1_r7.all_models.gff3.gz \
      osa1_r7.locus_brief_info.txt.gz NC_001320.1.fa NC_011033.1.fa genome.fa ) \
    > "$REF/reference_sha256.txt"
cut -f1,2 "$REF/genome.fa.fai"

# name the scripts read, file, md5 of the published file
ESM="https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41477-024-01754-4/MediaObjects"
while read -r name file md5; do
    gi_fetch "$ESM/$file" "$SUPP/$name"
    [ "$(md5sum < "$SUPP/$name" | cut -c1-32)" = "$md5" ] \
        || { echo "!!! $name does not match the published file; delete it and run again"; exit 1; }
done <<'EOF'
tonosaki2024_SupplementaryData1_NGS.xlsx 41477_2024_1754_MOESM3_ESM.xlsx 3f2ab1593d9e24fee7f1057ad13043e6
tonosaki2024_SupplementaryData3_imprint_classes.xlsx 41477_2024_1754_MOESM5_ESM.xlsx 9a7cc93ba6b6803900c656aeafd3d378
EOF
