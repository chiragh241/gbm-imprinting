#!/usr/bin/env bash
# Arabidopsis datasets used as deposited (no reads are mapped):
#
#   Allele-resolved endosperm CG methylation, reciprocal Col-0 x C24, 7 DAP (GEO GSE280598).
#   The Col_spiked files are on TAIR10 coordinates; the Col-0 allele is maternal in Col x C24
#   and paternal in C24 x Col:
#     GSM8601909-11  <GSM>_ColxC24_<rep>_Col_spiked_CG.bed.gz  -> mat_r1..3.CG.bed.gz
#     GSM8601912-14  <GSM>_C24xCol_<rep>_Col_spiked_CG.bed.gz  -> pat_r1..3.CG.bed.gz
#   Wild-type Col-0 endosperm (7 DAP) and young leaf, same series, EM-seq, kept under their GEO
#   names (Supplemental Fig. S7a-d):
#     GSM8601882-84  <GSM>_wt_<rep>_Col_spiked_CG.bed.gz          -> gse280598_tissues/
#     GSM8601894-96  <GSM>_wt_<rep>_leaf_Col_spiked_CG.bed.gz     -> gse280598_tissues/
#   Allele-resolved endosperm CG methylation in dme and fie crosses (Ibarra et al. 2012, GSE38935):
#     GSM952440_TFH605.all.cg-{dmecol,ler}.gff.gz      dme Col x Ler
#     GSM952441_TFH44.all.cg-{col,ler}.gff.gz          WT Col x Ler
#     GSM952442_TFH113.all.cg-{col,fieler}.gff.gz      fie Ler x Col
#     GSM952443_TFH52.all.cg-{col,ler}.merged.gff.gz   WT Ler x Col
#   Published supplementary tables, saved under the names the scripts read:
#     pignatta2014_endo.xlsx       Pignatta et al. 2014, eLife 3:e03198, Figure 1 source data 2
#                                  (elife-03198-fig1-data2-v4.xlsx)
#     gehring2011_tableS6.xls      Gehring et al. 2011, PLoS One 6:e23687, Table S6
#                                  (pone.0023687.s006.xls)
#     wolff2011_tableS4_MEG.xlsx   Wolff et al. 2011, PLoS Genetics 7:e1002126, Table S4
#                                  (pgen.1002126.s013.xlsx)
#     wolff2011_tableS7_PEG.xlsx   Wolff et al. 2011, Table S7 (pgen.1002126.s016.xlsx)
#     pnas.1604666113.sd01.xlsx    Bewick et al. 2016, PNAS 113:9111, Dataset S1. PNAS and PMC do
#                                  not serve it to scripts, so it is downloaded by hand (below).
#
# Output  $GI_DATA/arabidopsis/{gse280598,gse280598_tissues,gse38935,supplementary}/...
set -euo pipefail
source "$(dirname "$0")/../../config.sh"
A=$GI_DATA/arabidopsis
S=$A/supplementary
mkdir -p "$A/gse280598" "$A/gse280598_tissues" "$A/gse38935" "$S"
GEO=https://ftp.ncbi.nlm.nih.gov/geo/samples

BEW=$S/pnas.1604666113.sd01.xlsx
[ -s "$BEW" ] || {
    echo "!!! Download Dataset S1 of Bewick et al. 2016 in a browser (Supporting Information of"
    echo "    https://doi.org/10.1073/pnas.1604666113, file pnas.1604666113.sd01.xlsx),"
    echo "    save it as $BEW and run this script again."
    exit 1
}
[ "$(md5sum < "$BEW" | cut -c1-32)" = b5d7180b4ab80fa95d57619ba4c96802 ] \
    || { echo "!!! $BEW does not match the published Dataset S1"; exit 1; }

# name the scripts read, URL, md5 of the published file
while read -r name url md5; do
    gi_fetch "$url" "$S/$name"
    [ "$(md5sum < "$S/$name" | cut -c1-32)" = "$md5" ] \
        || { echo "!!! $name does not match the published file; delete it and run again"; exit 1; }
done <<'EOF'
pignatta2014_endo.xlsx https://cdn.elifesciences.org/articles/03198/elife-03198-fig1-data2-v4.xlsx 23b73f67b673ab28ff4024f3607d26aa
gehring2011_tableS6.xls https://journals.plos.org/plosone/article/file?type=supplementary&id=10.1371/journal.pone.0023687.s006 f10c4f73a63ce86b285cbeb5823e3ba9
wolff2011_tableS4_MEG.xlsx https://journals.plos.org/plosgenetics/article/file?type=supplementary&id=10.1371/journal.pgen.1002126.s013 9b71e0e33d5d216b5c925bd23e5d1eab
wolff2011_tableS7_PEG.xlsx https://journals.plos.org/plosgenetics/article/file?type=supplementary&id=10.1371/journal.pgen.1002126.s016 77e85fb09f677c5acfc04f6806ac6e2b
EOF

geo_file() {  # GSM, file-name regex, output path
    local gsm=$1 pat=$2 out=$3 dir fn
    [ -s "$out" ] && return 0
    dir=$GEO/${gsm%???}nnn/$gsm/suppl
    fn=$(curl -s "$dir/" | grep -oE "$pat" | sort -u | head -1)
    [ -n "$fn" ] || { echo "!!! nothing matching $pat for $gsm"; return 1; }
    gi_fetch "$dir/$fn" "$out"
    echo "$gsm: $fn -> $(basename "$out")"
}

rep=1
for gsm in GSM8601909 GSM8601910 GSM8601911; do
    geo_file "$gsm" "${gsm}_ColxC24_[0-9]+_Col_spiked_CG\.bed\.gz" "$A/gse280598/mat_r${rep}.CG.bed.gz"
    rep=$((rep + 1))
done
rep=1
for gsm in GSM8601912 GSM8601913 GSM8601914; do
    geo_file "$gsm" "${gsm}_C24xCol_[0-9]+_Col_spiked_CG\.bed\.gz" "$A/gse280598/pat_r${rep}.CG.bed.gz"
    rep=$((rep + 1))
done

for f in GSM952440_TFH605.all.cg-dmecol.gff.gz GSM952440_TFH605.all.cg-ler.gff.gz \
         GSM952441_TFH44.all.cg-col.gff.gz GSM952441_TFH44.all.cg-ler.gff.gz \
         GSM952442_TFH113.all.cg-col.gff.gz GSM952442_TFH113.all.cg-fieler.gff.gz \
         GSM952443_TFH52.all.cg-col.merged.gff.gz GSM952443_TFH52.all.cg-ler.merged.gff.gz; do
    gsm=${f%%_*}
    gi_fetch "$GEO/${gsm%???}nnn/$gsm/suppl/$f" "$A/gse38935/$f"
done

for f in GSM8601882_wt_1_Col_spiked_CG.bed.gz GSM8601883_wt_2_Col_spiked_CG.bed.gz \
         GSM8601884_wt_3_Col_spiked_CG.bed.gz GSM8601894_wt_1_leaf_Col_spiked_CG.bed.gz \
         GSM8601895_wt_2_leaf_Col_spiked_CG.bed.gz GSM8601896_wt_3_leaf_Col_spiked_CG.bed.gz; do
    gsm=${f%%_*}
    gi_fetch "$GEO/${gsm%???}nnn/$gsm/suppl/$f" "$A/gse280598_tissues/$f"
done

# GEO publishes no checksums for supplementary files; the gzip check catches a truncated one
for f in "$A"/gse280598/*.gz "$A"/gse280598_tissues/*.gz "$A"/gse38935/*.gz; do gzip -t "$f" || { echo "!!! corrupt: $f"; exit 1; }; done
ls -la "$A" "$A/gse280598" "$A/gse280598_tissues" "$A/gse38935"
