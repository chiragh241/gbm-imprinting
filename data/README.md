# Data files shipped with the repository

Small inputs that are not reliably downloadable by script, or are derived tables whose
derivation is documented here. Everything else is downloaded by the numbered scripts (see the
top-level README), and the large processed files are on Zenodo.

## maize/

| File | What it is | Source |
|---|---|---|
| `B73_NAM_syntelog_maternalpref.txt.gz` | Allele-specific endosperm expression (RPM) of B73 reciprocal hybrids with eight NAM founders at 14 DAP, one row per gene, hybrid and alignment reference, with each row's B73 v5 syntelog(s) | 237,338 of the 444,878 rows of `mat_pref_NAM.txt.gz`, Higgins et al. 2025 (G3), https://github.com/kmhiggins/Imprinting_timeseries_2024, allele values unchanged, with two added columns: `feature_2`, the Zm00001eb genes in the gene's row of the MaizeGDB 2020 pan-gene table (for B73-aligned rows, the gene itself), and `cont_genome` (`<genome> _ <contrast>`). The kept rows are those with A_mean + B_mean of at least 2 RPM, B73-syntelog rows only for founder-aligned rows; that rule reproduces all 237,338 rows but would add 848 more, so the table is shipped as it was used. |
| `B73_endo-leaf_hyperDMRs.bed` | 27,402 B73 v5 regions of 200 bp less CG-methylated in endosperm than in leaf. Columns: 1-3 the region (chromosome, start, end); 4 CG methylation in endosperm; 5 CG methylation in leaf; 6 read coverage in endosperm; 7 read coverage in leaf; 8 number of informative cytosines (covered by reads in endosperm). The analysis uses columns 1-3. | Provided by Jonathan Gent; not published elsewhere. |

## emseq/, allele/, rnaseq/

Run lists: SRA run accessions with sample tags (`emseq/runs.tsv`, `allele/runs.tsv`), and the
sample sheet and ENA FASTQ md5 checksums of the Xu et al. 2022 RNA-seq (`rnaseq/`).

## arabidopsis

Nothing is shipped. The published supplementary tables (Pignatta et al. 2014, Gehring et al. 2011, Wolff et al. 2011,
Bewick et al. 2016) are downloaded by
`scripts/11_arabidopsis/01_download.sh`, which lists each one's source file and md5. Bewick et
al. 2016 Dataset S1 has to be downloaded by hand in a browser; the script says where to save it.

## rice/

| File | What it is | Source |
|---|---|---|
| `runs.tsv` | The eight 5 DAP F1 endosperm WGBS libraries of Tonosaki et al. 2024 (DRA015356) with their cross and replicate, and the four Kasalath and 93-11 resequencing runs they list (Supplementary Data 1); cross directions confirmed from the BioSample cultivar and library names | Compiled for this work |
| `rice_leaf_methylation.csv` | Per-gene CG, CHG and CHH methylation of Nipponbare young leaf (SRR10763664, Zhao et al. 2020, Nat Commun 11:2658) for 38,673 MSU7 genes (38,672 of the 39,049 not annotated as transposable element related, and one that is): methylated, unmethylated and total read counts over each context's sites (`*_mC`, `*_umC`, `*_totalC`), their ratios (`*_methy`), and the category (`epiallele`), which `13_rice/06_gene_sets.py` checks against the rules it applies. sha256 9890849e424f1a03421b2cc7e3cd55cc824474531eaacf84c7bb3e69ada32862 | Provided by Yirui Sun; not published elsewhere |

Tonosaki et al.'s Supplementary Data 1 and 3 are downloaded by `scripts/13_rice/01_download_references.sh`,
which checks each one's md5.
