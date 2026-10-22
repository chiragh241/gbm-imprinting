#!/usr/bin/env python
"""Per-library QC of the eight rice WGBS libraries beside Tonosaki et al.'s own numbers
(Supplementary Data 1, WGBS sheet), from the reports 05_map_libraries.sh keeps.

  raw pairs          cutadapt "Total read pairs processed"        SD1 "Raw reads pair"
  unique pairs       Bismark "unique best hit"                    SD1 "Uniquly mapped reads pair"
  dedup pairs        left after deduplication (all chromosomes)
  Nip / alt pairs    SNPsplit sort report, read pairs per genome
  Nip / alt reads    SNPsplit tag report, single reads per genome  SD1 "Allelic reads pair",
                     (SD1's allele columns match these magnitudes,   reference / alternative
                     though its header says pairs)
  mat / pat          maternal over paternal allele pairs (2 for unbiased endosperm)
  conversion         chloroplast reads: 100 x (1 - methylated / all C calls), all contexts and
                     CHH alone                                     SD1 "Conversion rate"

Input   $GI_WORK/13_rice/qc/   (05)
        $GI_DATA/rice/supplementary/tonosaki2024_SupplementaryData1_NGS.xlsx   (01)
Output  $GI_RESULTS/tables/rice_qc_vs_tonosaki.csv
"""
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, WORK, SHIPPED, TAB

Q = WORK / "13_rice" / "qc"


def grab(path, pattern):
    m = re.search(pattern, path.read_text())
    if not m:
        raise SystemExit(f"{pattern!r} not in {path.name}")
    return int(m.group(1).replace(",", ""))


runs = pd.read_csv(SHIPPED / "rice" / "runs.tsv", sep="\t")
runs = runs[runs.type == "wgbs"]
sd1 = pd.read_excel(DATA / "rice" / "supplementary" / "tonosaki2024_SupplementaryData1_NGS.xlsx",
                    sheet_name="WGBS", header=None, skiprows=4).iloc[:, [5, 6, 8, 9, 10, 11]]
sd1.columns = ["run", "sd1_raw_pairs", "sd1_unique_pairs", "sd1_conversion", "sd1_nip", "sd1_alt"]
sd1 = sd1[sd1.run.isin(runs.run)].set_index("run")

rows = []
for r in runs.itertuples():
    stem = f"{r.run}_1_val_1_bismark_bt2_pe"
    chrc = Q / f"{r.run}.chrC_splitting_report.txt"
    meth = {c: grab(chrc, rf"Total methylated C's in {c} context:\s+(\d+)") for c in ("CpG", "CHG", "CHH")}
    conv = {c: grab(chrc, rf"Total C to T conversions in {c} context:\s+(\d+)") for c in ("CpG", "CHG", "CHH")}
    nip = grab(Q / f"{stem}.deduplicated.SNPsplit_sort.txt", r"specific for genome 1:\s+(\d+)")
    alt = grab(Q / f"{stem}.deduplicated.SNPsplit_sort.txt", r"specific for genome 2:\s+(\d+)")
    rows.append(dict(
        run=r.run, mother=r.mother, father=r.father, rep=r.rep,
        raw_pairs=grab(Q / f"{r.run}_cutadapt_report.txt", r"Total read pairs processed:\s+([\d,]+)"),
        unique_pairs=grab(Q / f"{r.run}_1_val_1_bismark_bt2_PE_report.txt", r"unique best hit:\s+(\d+)"),
        dedup_pairs=grab(Q / f"{stem}.deduplication_report.txt", r"duplicates removed \([\d.]+%\), (\d+) left"),
        nip_pairs=nip, alt_pairs=alt,
        conflicting_pairs=grab(Q / f"{stem}.deduplicated.SNPsplit_sort.txt",
                               r"conflicting SNP information:\s+(\d+)"),
        nip_reads=grab(Q / f"{stem}.deduplicated.SNPsplit_report.txt", r"(\d+) reads were specific for genome 1"),
        alt_reads=grab(Q / f"{stem}.deduplicated.SNPsplit_report.txt", r"(\d+) reads were specific for genome 2"),
        mat_over_pat=(nip / alt) if r.mother == "NIP" else (alt / nip),
        conversion_all=100 * (1 - sum(meth.values()) / (sum(meth.values()) + sum(conv.values()))),
        conversion_chh=100 * (1 - meth["CHH"] / (meth["CHH"] + conv["CHH"])),
    ))
t = pd.DataFrame(rows).set_index("run").join(sd1)
t["unique_over_sd1"] = t.unique_pairs / t.sd1_unique_pairs
t["nip_reads_over_sd1"] = t.nip_reads / t.sd1_nip
t["alt_reads_over_sd1"] = t.alt_reads / t.sd1_alt
TAB.mkdir(parents=True, exist_ok=True)
t.to_csv(TAB / "rice_qc_vs_tonosaki.csv", float_format="%.4f")
with pd.option_context("display.width", 250, "display.max_columns", 30):
    print(t.round(3).to_string())
