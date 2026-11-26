"""Per-gene replication-timing (RT) score from one smoothed hg19 RT bigWig.

Used for every track: the ENCODE wavelet-smoothed tracks and the drug-treatment
tracks of wavelet_rt.py. For every protein-coding gene:
  1. Gene signal. Mean value of the track intervals that overlap the gene body
     [start, end]. If no interval overlaps the gene (a short gene in a sparse
     track), the signal is interpolated linearly at the gene-body midpoint from
     the nearest interval on each side, if these are at most 100 kb apart.
  2. Orientation. If the mean signal of 11 early housekeeping genes is lower than
     that of 6 late anchor genes, the signal is sign-flipped (high = early).
  3. Genome reference. Mean value of the intervals in each consecutive 10 kb bin
     over chr1-22 + chrX.
  4. rt_score = 1 - mid-rank percentile of the gene signal among the genome
     bins. 0 = earliest, 1 = latest.

Usage: gene_rt.py <track.bw> <genes.parquet> <out.parquet> [track_id]
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
BIN = 10_000
MAX_GAP = 100_000   # widest data gap bridged by interpolation


def read_track(path):
    """{chrom: (starts, ends, values, length)} of the finite intervals of a bigWig."""
    bw = pyBigWig.open(str(path))
    out = {}
    for c, clen in bw.chroms().items():
        iv = np.array(bw.intervals(c) or [], dtype=float).reshape(-1, 3)
        iv = iv[np.isfinite(iv[:, 2])]
        out[c] = (iv[:, 0], iv[:, 1], iv[:, 2], int(clen))
    bw.close()
    return out


def gene_signal(track, chrom, start, end):
    sig = np.full(len(chrom), np.nan)
    for i, (c, s, e) in enumerate(zip(chrom, start, end)):
        if c not in track or len(track[c][0]) == 0:
            continue
        st, en, v, _ = track[c]
        lo = np.searchsorted(en, s, "right")   # first interval ending after the gene start
        hi = np.searchsorted(st, e, "left")    # first interval starting at or after the gene end
        if hi > lo:
            sig[i] = v[lo:hi].mean()
        elif 0 < lo < len(st) and st[lo] - en[lo - 1] <= MAX_GAP:   # no overlap: interpolate
            mid = (s + e) / 2
            x0, x1 = (st[lo - 1] + en[lo - 1]) / 2, (st[lo] + en[lo]) / 2
            sig[i] = v[lo - 1] + (v[lo] - v[lo - 1]) * (mid - x0) / (x1 - x0)
    return sig


def genome_bins(track):
    out = []
    for c in CHROMS:
        if c not in track:
            continue
        st, _, v, clen = track[c]
        nb = clen // BIN
        b = (st // BIN).astype(int)
        keep = b < nb
        n = np.bincount(b[keep], minlength=nb)
        tot = np.bincount(b[keep], weights=v[keep], minlength=nb)
        out.append(tot[n > 0] / n[n > 0])
    return np.concatenate(out)


def to_rt(sig, g, is_hk, is_late):
    flip = -1.0 if np.nanmean(sig[is_hk]) < np.nanmean(sig[is_late]) else 1.0
    gs = np.sort(g * flip)
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
    ap.add_argument("track_id", nargs="?")
    a = ap.parse_args()

    genes = pd.read_parquet(a.genes).sort_values(["chrom", "start"]).reset_index(drop=True)
    track = read_track(a.track)
    g = genome_bins(track)
    sig = gene_signal(track, genes.chrom.values, genes.start.values, genes.end.values)
    out = genes[["gene_id", "gene_name"]].copy()
    out["signal"] = sig
    out["rt_score"], out["flip"] = to_rt(sig, g, genes.gene_name.isin(HOUSEKEEPING).values,
                                         genes.gene_name.isin(LATE_REF).values)
    out["track"] = a.track_id or Path(a.track).name.split(".")[0]
    out.to_parquet(a.out, index=False)
    print(f"{a.track}: {out.rt_score.notna().sum()} genes with a value; {len(g)} genome bins; flip={out.flip.iat[0]}")


if __name__ == "__main__":
    main()
