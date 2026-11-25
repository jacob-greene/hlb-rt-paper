"""Flag the top 5% expressed protein-coding genes per cell line from ENCODE RNA-seq.

Per cell line: TPM is averaged over all its gene-quantification files (one per
replicate), per unversioned Ensembl gene ID. The table is restricted to the
protein-coding genes of all_coding.parquet. A gene is "top 5%" if its mean TPM
is at or above the 95th percentile of that cell line's protein-coding genes.

Usage: rnaseq_top5.py <samples_rnaseq.csv> <tsv_dir> <all_coding.parquet> <out.parquet>
"""
import sys
from pathlib import Path

import pandas as pd

sheet, tsv_dir, coding, out = sys.argv[1:5]
pc_ids = set(pd.read_parquet(coding).gene_id_unversioned)
s = pd.read_csv(sheet)
# The pipeline can run on a subset of cell lines: keep the files that were fetched.
s = s[[(Path(tsv_dir) / f"{a}.tsv").exists() for a in s.file_accession]]
assert len(s), "no RNA-seq file of the sheet is in tsv_dir"
rows = []
for cell, sub in s.groupby("cell_line", sort=False):
    parts = []
    for acc in sub.file_accession:
        df = pd.read_csv(Path(tsv_dir) / f"{acc}.tsv", sep="\t", low_memory=False)
        cols = {c.lower(): c for c in df.columns}
        tpm = cols.get("tpm") or cols.get("expected_tpm")
        parts.append(df[[cols["gene_id"], tpm]].set_axis(["gene_id", "tpm"], axis=1))
    e = pd.concat(parts, ignore_index=True)
    e["gene_id_unversioned"] = e.gene_id.str.split(".").str[0]
    avg = e.groupby("gene_id_unversioned").tpm.mean().reset_index()
    avg = avg[avg.gene_id_unversioned.isin(pc_ids)].copy()
    avg = avg[avg.tpm >= avg.tpm.quantile(0.95)]
    rows.append(avg.assign(cell_line=cell))
top5 = pd.concat(rows, ignore_index=True)[["cell_line", "gene_id_unversioned", "tpm"]]
top5.to_parquet(out, index=False)
print(f"top-5% rows: {len(top5)}; cell lines: {top5.cell_line.nunique()}")
