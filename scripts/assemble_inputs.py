"""Collect the pipeline outputs into the three tables that analysis.ipynb reads.

Layout of <outdir>/ (the same relative paths as the original analysis):
  aggregated/per_file_rtscore.parquet  one row per (track x gene) with an RT score
  genesets/histones.parquet            replication-dependent histone genes
  rnaseq/top5_per_cell.parquet         top-5% expressed protein-coding genes per label

Usage: assemble_inputs.py <sheet.csv> <histones.parquet> <top5.parquet> <outdir> <gene_rt.parquet> ...
"""
import shutil
import sys
from pathlib import Path

import pandas as pd

sheet, hist, top5, outdir, *files = sys.argv[1:]
out = Path(outdir)
for sub in ("aggregated", "genesets", "rnaseq"):
    (out / sub).mkdir(parents=True, exist_ok=True)

s = pd.read_csv(sheet).drop_duplicates("file_accession").set_index("file_accession")
rt = pd.concat([pd.read_parquet(f, columns=["track", "gene_id", "rt_tss5kb"]) for f in files], ignore_index=True)
rt = rt.rename(columns={"track": "file_accession", "rt_tss5kb": "rt_score"}).dropna(subset=["rt_score"])
for col in ("cell_line", "group", "assay"):
    rt[col] = rt.file_accession.map(s[col])
assert rt.cell_line.notna().all(), "track missing from sample sheet"
rt = rt[["file_accession", "cell_line", "group", "assay", "gene_id", "rt_score"]]
rt.to_parquet(out / "aggregated/per_file_rtscore.parquet", index=False)
shutil.copy(hist, out / "genesets/histones.parquet")
shutil.copy(top5, out / "rnaseq/top5_per_cell.parquet")
print(f"{rt.file_accession.nunique()} tracks, {rt.cell_line.nunique()} labels, {len(rt)} rows -> {out}")
