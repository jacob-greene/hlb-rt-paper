"""Per-gene replication-timing (RT) score from one smoothed hg19 RT bigWig (ENCODE wavelet track).

For every protein-coding gene: the mean signal over TSS +/- 5 kb, converted to an RT score:
  1. Orientation. If the mean of 11 early housekeeping genes is lower than the
     mean of 6 late anchor genes, the signal is sign-flipped (high = early).
  2. Genome reference. Mean signal in consecutive 10 kb bins over chr1-22 + chrX.
  3. rt_score = 1 - mid-rank percentile of the gene signal among the genome
     bins. 0 = earliest, 1 = latest.

Signal means use pyBigWig ``stats(type="mean")`` (zoom-level summaries where
available), exactly as in the original analysis.

Usage: gene_rt.py <track.bw> <genes.parquet> <out.parquet>
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
HOUSEKEEPING = ["ACTB", "GAPDH", "MYC", "POLR2A", "EEF1A1", "RPS18",
                "RPL18", "TBP", "PPIA", "B2M", "HPRT1"]
LATE_REF = ["DLG2", "PCDH7", "OLFM3", "ROBO2", "CSMD1", "AFM"]
TSS_WIN = 5000
BIN = 10_000


def window_means(bw, chrom, start, end):
    ch = bw.chroms()
    sig = np.full(len(chrom), np.nan)
    for i, (c, s, e) in enumerate(zip(chrom, start, end)):
        if c not in ch or not np.isfinite(s):
            continue
        s, e = max(0, int(s)), min(int(e), int(ch[c]))
        if e <= s:
            continue
        v = bw.stats(c, s, e, type="mean")[0]
        sig[i] = np.nan if v is None else float(v)
    return sig


def genome_bins(bw):
    ch = bw.chroms()
    out = []
    for c in CHROMS:
        if c in ch:
            nb = max(1, ch[c] // BIN)
            out.append(np.array([np.nan if x is None else float(x)
                                 for x in bw.stats(c, 0, ch[c], nBins=nb, type="mean")]))
    return np.concatenate(out)


def to_rt(sig, g, is_hk, is_late):
    flip = -1.0 if np.nanmean(sig[is_hk]) < np.nanmean(sig[is_late]) else 1.0
    gs = np.sort(g[~np.isnan(g)] * flip)
    x = sig * flip
    rt = np.full(x.shape, np.nan)
    ok = ~np.isnan(x)
    rt[ok] = 1 - (np.searchsorted(gs, x[ok], "left") + np.searchsorted(gs, x[ok], "right")) / 2 / len(gs)
    return rt, flip


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("track")
    ap.add_argument("genes")
    ap.add_argument("out")
    a = ap.parse_args()

    genes = pd.read_parquet(a.genes)
    genes = genes[["gene_id", "gene_name", "chrom", "start", "end", "strand", "tss"]].copy()

    bw = pyBigWig.open(a.track)
    g = genome_bins(bw)
    is_hk = genes.gene_name.isin(HOUSEKEEPING).values
    is_late = genes.gene_name.isin(LATE_REF).values
    out = genes[["gene_id", "gene_name"]].copy()
    s, e = (genes.tss - TSS_WIN).clip(lower=0), genes.tss + TSS_WIN
    sig = window_means(bw, genes.chrom.values, s.values, e.values)
    out["rt_tss5kb"], out["flip"] = to_rt(sig, g, is_hk, is_late)
    bw.close()
    out["track"] = Path(a.track).name.split(".")[0]
    out.to_parquet(a.out, index=False)
    print(f"{a.track}: {out.rt_tss5kb.notna().sum()} genes with a value; {len(g)} genome bins; flip={out.flip.iat[0]}")


if __name__ == "__main__":
    main()
