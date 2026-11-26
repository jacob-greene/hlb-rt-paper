"""Collect the per-gene RT scores (gene_rt.py) into the tables that analysis.ipynb reads.

Layout of <outdir>/ (the same relative paths as the original analysis):
  aggregated/per_file_rtscore.parquet  one row per (ENCODE track x gene) with an RT score
  genesets/histones.parquet            replication-dependent histone genes
  rnaseq/top5_per_cell.parquet         top-5% expressed protein-coding genes per label
  drug/drug_per_gene.csv               one row per (drug condition x histone gene) with an RT score

Usage: assemble_inputs.py <sheet.csv> <histones.parquet> <top5.parquet> <outdir>
           --encode <gene_rt.parquet> ... --drug <condition>=<gene_rt.parquet> ...
"""
import argparse
import shutil
from pathlib import Path

import pandas as pd


def assign_class(sym):
    s = str(sym).upper()
    if s.startswith("HIST"):
        rest = s.split("H", 2)
        if len(rest) >= 3:
            head = rest[2]
            if head.startswith("1"): return "H1"
            if head.startswith("2A"): return "H2A"
            if head.startswith("2B"): return "H2B"
            if head.startswith("3"): return "H3"
            if head.startswith("4"): return "H4"
    if s.startswith("H1-"): return "H1"
    if s.startswith("H2AC") or s.startswith("H2A"): return "H2A"
    if s.startswith("H2BC"): return "H2B"
    if s.startswith("H3C") or s.startswith("H3-"): return "H3"
    if s.startswith("H4C") or s.startswith("H4-"): return "H4"
    return "variant"


ap = argparse.ArgumentParser()
ap.add_argument("sheet")
ap.add_argument("histones")
ap.add_argument("top5")
ap.add_argument("outdir")
ap.add_argument("--encode", nargs="+", required=True)
ap.add_argument("--drug", nargs="+", required=True)
a = ap.parse_args()
out = Path(a.outdir)
for sub in ("aggregated", "genesets", "rnaseq", "drug"):
    (out / sub).mkdir(parents=True, exist_ok=True)

s = pd.read_csv(a.sheet).drop_duplicates("file_accession").set_index("file_accession")
rt = pd.concat([pd.read_parquet(f, columns=["track", "gene_id", "rt_score"]) for f in a.encode], ignore_index=True)
rt = rt.rename(columns={"track": "file_accession"}).dropna(subset=["rt_score"])
for col in ("cell_line", "group", "assay"):
    rt[col] = rt.file_accession.map(s[col])
assert rt.cell_line.notna().all(), "track missing from sample sheet"
rt = rt[["file_accession", "cell_line", "group", "assay", "gene_id", "rt_score"]]
rt.to_parquet(out / "aggregated/per_file_rtscore.parquet", index=False)

hist = pd.read_parquet(a.histones)
drug = []
for spec in a.drug:                      # in condition order
    condition, path = spec.split("=", 1)
    d = pd.read_parquet(path).merge(hist[["gene_id", "histone_cluster"]], on="gene_id")
    d = d.dropna(subset=["rt_score"]).assign(condition=condition)
    d["histone_class"] = d.gene_name.map(assign_class)
    drug.append(d.rename(columns={"signal": "raw_signal"})[
        ["gene_id", "gene_name", "histone_cluster", "histone_class", "condition", "raw_signal", "rt_score"]])
drug = pd.concat(drug, ignore_index=True)
drug.to_csv(out / "drug/drug_per_gene.csv", index=False)

shutil.copy(a.histones, out / "genesets/histones.parquet")
shutil.copy(a.top5, out / "rnaseq/top5_per_cell.parquet")
print(f"{rt.file_accession.nunique()} tracks, {rt.cell_line.nunique()} labels, {len(rt)} rows; "
      f"drug: {len(drug)} rows, {drug.condition.nunique()} conditions -> {out}")
