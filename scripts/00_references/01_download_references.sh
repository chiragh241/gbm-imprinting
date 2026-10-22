#!/usr/bin/env bash
# Genome assemblies, gene annotations and gene cross-references used by several modules.
#
# Every file is public; each is fetched once into $GI_DATA/references and never modified in place.
#
#   B73 v5 genome        Ensembl Plants release 57
#                        https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/zea_mays/dna/
#                        Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz   -> B73v5.fa
#                        RNA-seq and ChIP-seq alignment (chromosomes named 1..10, no "chr")
#   B73 v5 genes         MaizeGDB, https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/
#                        Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz            -> B73v5.gff3,
#                        and its gene lines alone -> B73v5_genes.gff3 (gene coordinates, "chr1..")
#   W22 v2 genome, genes MaizeGDB, https://download.maizegdb.org/Zm-W22-REFERENCE-NRGENE-2.0/
#                        Zm-W22-REFERENCE-NRGENE-2.0.fa.gz                       -> W22v2.fa
#                        Zm-W22-REFERENCE-NRGENE-2.0_Zm00004b.1.gff3.gz          -> W22v2.gff3.gz
#   Mo17 CAU-2.0 genome  MaizeGDB, https://download.maizegdb.org/Zm-Mo17-REFERENCE-CAU-2.0/
#                        Zm-Mo17-REFERENCE-CAU-2.0.fa.gz                         -> Mo17_CAU2.fa.gz
#   B73 v4 genes         MaizeGDB, https://download.maizegdb.org/Zm-B73-REFERENCE-GRAMENE-4.0/
#                        Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3.gz         -> B73v4.gff3.gz
#   B73 v5 to v4 genes   MaizeGDB, https://download.maizegdb.org/Pan-genes/B73_gene_xref/
#                        B73v5_to_B73v4.tsv
#   Maize pan-genes      MaizeGDB, https://download.maizegdb.org/Pan-genes/archive/pan-zea-2020/
#                        MaizeGDB_maize_pangene_2020_08.tsv.gz  (B73 v5 to W22 v2 syntelogs)
#   Arabidopsis genes    Ensembl Plants release 57
#                        https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/gff3/arabidopsis_thaliana/
#                        Arabidopsis_thaliana.TAIR10.57.gff3.gz
set -euo pipefail
source "$(dirname "$0")/../../config.sh"
R=$GI_DATA/references
mkdir -p "$R"

ENS=https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57
MGDB=https://download.maizegdb.org

gi_fetch "$ENS/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz" \
    "$R/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz"
gi_fetch "$MGDB/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz" \
    "$R/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz"
gi_fetch "$MGDB/Zm-W22-REFERENCE-NRGENE-2.0/Zm-W22-REFERENCE-NRGENE-2.0.fa.gz" \
    "$R/Zm-W22-REFERENCE-NRGENE-2.0.fa.gz"
gi_fetch "$MGDB/Zm-W22-REFERENCE-NRGENE-2.0/Zm-W22-REFERENCE-NRGENE-2.0_Zm00004b.1.gff3.gz" \
    "$R/W22v2.gff3.gz"
gi_fetch "$MGDB/Zm-Mo17-REFERENCE-CAU-2.0/Zm-Mo17-REFERENCE-CAU-2.0.fa.gz" \
    "$R/Mo17_CAU2.fa.gz"
gi_fetch "$MGDB/Zm-B73-REFERENCE-GRAMENE-4.0/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3.gz" \
    "$R/B73v4.gff3.gz"
gi_fetch "$MGDB/Pan-genes/B73_gene_xref/B73v5_to_B73v4.tsv" "$R/B73v5_to_B73v4.tsv"
gi_fetch "$MGDB/Pan-genes/archive/pan-zea-2020/MaizeGDB_maize_pangene_2020_08.tsv.gz" \
    "$R/MaizeGDB_maize_pangene_2020_08.tsv.gz"
gi_fetch "$ENS/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.57.gff3.gz" \
    "$R/Arabidopsis_thaliana.TAIR10.57.gff3.gz"

for f in "$R"/*.gz; do gzip -t "$f" || { echo "corrupt: $f" >&2; exit 1; }; done

# unpacked under a temporary name, so an interrupted run never leaves a truncated file to skip
unpack() { [ -s "$2" ] || { zcat "$1" > "$2.part"; mv -f "$2.part" "$2"; }; }
unpack "$R/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz" "$R/B73v5.fa"
unpack "$R/Zm-W22-REFERENCE-NRGENE-2.0.fa.gz" "$R/W22v2.fa"
unpack "$R/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz" "$R/B73v5.gff3"
awk -F'\t' '$3 == "gene"' "$R/B73v5.gff3" > "$R/B73v5_genes.gff3"

echo "B73 v5 contigs: $(grep -c '^>' "$R/B73v5.fa")   W22 v2 contigs: $(grep -c '^>' "$R/W22v2.fa")"
echo "B73 v5 genes: $(wc -l < "$R/B73v5_genes.gff3")"
