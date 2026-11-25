"""Per-gene replication-timing (RT) score of the histone genes in the drug-treatment tracks.

Port of 16b_drug_wavelet_exact.py of the original analysis (aggregation and table only).
Input tracks: the output of wavelet_rt.py for each condition (1 kb bins, high = early).

For every protein-coding gene:
  1. signal = mean of the 1 kb bins that overlap the exact gene body [start, end];
  2. genome reference = mean signal in consecutive 10 kb bins over chr1-22 + chrX;
  3. rt_score = 1 - mid-rank percentile of the gene signal among the genome bins.
     0 = earliest, 1 = latest.
Only the histone genes are written out.

Usage: drug_gene_rt.py <all_coding.parquet> <histones.parquet> <out.csv> <condition>=<track.bw> ...
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
BIN_BP = 1000


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


def aggregate_one(bw_path, chrom, start, end):
    """Gene-body mean signal and genome-rank RT score for each interval."""
    bw = pyBigWig.open(str(bw_path))
    lens = bw.chroms()
    arrs = {}
    for c, clen in lens.items():
        nb = int(clen) // BIN_BP
        if nb == 0:
            continue
        a = np.full(nb, np.nan)
        for s, e, v in bw.intervals(c) or []:
            i = s // BIN_BP
            if 0 <= i < nb:
                a[i] = v
        arrs[c] = a
    bw.close()

    sig = np.full(len(chrom), np.nan)
    for i, (c, s, e) in enumerate(zip(chrom, start, end)):
        if c not in arrs:
            continue
        a = arrs[c]
        s_bp, e_bp = max(0, int(s)), min(int(e), len(a) * BIN_BP)
        if e_bp <= s_bp:
            continue
        s_bin, e_bin = s_bp // BIN_BP, (e_bp + BIN_BP - 1) // BIN_BP
        if s_bin >= len(a):
            continue
        slc = a[s_bin:min(e_bin, len(a))]
        slc = slc[np.isfinite(slc)]
        if len(slc):
            sig[i] = float(slc.mean())

    genome = []
    for c in CHROMS:
        if c not in arrs:
            continue
        a = arrs[c]
        nb10 = len(a) // 10
        if nb10 == 0:
            continue
        with warnings.catch_warnings():   # all-NaN 10 kb bins (gaps) give NaN, dropped below
            warnings.simplefilter("ignore", RuntimeWarning)
            genome.append(np.nanmean(a[:nb10 * 10].reshape(nb10, 10), axis=1))
    g = np.concatenate(genome) if genome else np.array([])
    g = np.sort(g[~np.isnan(g)])
    rt = np.full(len(chrom), np.nan)
    if len(g) < 100:
        return sig, rt
    ok = ~np.isnan(sig)
    lo = np.searchsorted(g, sig[ok], side="left")
    hi = np.searchsorted(g, sig[ok], side="right")
    rt[ok] = 1.0 - (lo + hi) / 2.0 / len(g)
    return sig, rt


def main():
    coding_path, histones_path, out_path, *tracks = sys.argv[1:]
    pc = pd.read_parquet(coding_path).sort_values(["chrom", "start"]).reset_index(drop=True)
    histones = pd.read_parquet(histones_path).set_index("gene_id")
    is_hist = pc["gene_id"].isin(histones.index).to_numpy()

    rows = []
    for spec in tracks:
        condition, path = spec.split("=", 1)
        sig, rt = aggregate_one(path, pc["chrom"].tolist(), pc["start"].tolist(), pc["end"].tolist())
        print(f"{condition}: {np.isfinite(sig).sum()} genes with signal, "
              f"all-coding median rt = {np.nanmedian(rt):.3f}")
        for i in np.flatnonzero(is_hist & np.isfinite(rt)):
            gid = pc.at[i, "gene_id"]
            h = histones.loc[gid]
            rows.append({
                "gene_id": gid,
                "gene_name": h["gene_name"],
                "histone_cluster": h["histone_cluster"],
                "histone_class": assign_class(h["gene_name"]),
                "condition": condition,
                "track": Path(path).name,
                "raw_signal": float(sig[i]),
                "rt_score": float(rt[i]),
            })
    out = pd.DataFrame(rows)
    out.to_csv(out_path, index=False)
    print(f"wrote {out_path}: {len(out)} rows, {out['gene_id'].nunique()} genes")


if __name__ == "__main__":
    main()
