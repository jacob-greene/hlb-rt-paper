"""Build one replication-timing track from six cell-cycle fraction bigWigs (Repli-seq, 1 kb bins).

Port of 15_wavelet_smooth.py of the original analysis:
  1. read G1b, S1, S2, S3, S4, G2 signal at 1 kb,
  2. percentage-normalise the six values at each bin,
  3. weighted average with the ENCODE/Thurman weights 0.917, 0.750, 0.583, 0.417, 0.250, 0,
  4. Daubechies-4 wavelet smoothing at level 7 (J7, about 128 kb), details set to zero.
Higher values = earlier replication. The output is on a 0-1 scale (the ENCODE track is 0-100).

The pipeline uses it twice:
  - drug-treatment tracks (Zenodo 10.5281/zenodo.19093559): input of the notebook, sections 3-5;
  - QC: rebuild of the 16 ENCODE UW tracks from their fraction signals (not used by the notebook).

Usage: wavelet_rt.py <out.bw> <G1b.bw> <S1.bw> <S2.bw> <S3.bw> <S4.bw> <G2.bw>
"""
import sys

import numpy as np
import pyBigWig
import pywt

WAVELET, LEVEL, BIN_BP = "db4", 7, 1000
WEIGHTS = np.array([0.917, 0.750, 0.583, 0.417, 0.250, 0.0])
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def wavelet_smooth_J7(sig):
    """db4 DWT at 7 levels; zero details; reconstruct. NaN bins are filled with the chromosome median."""
    valid = np.isfinite(sig)
    if valid.sum() == 0:
        return sig.copy()
    s = sig.copy()
    s[~valid] = float(np.median(sig[valid]))
    coeffs = pywt.wavedec(s, WAVELET, level=LEVEL, mode="symmetric")
    rec = pywt.waverec([coeffs[0]] + [np.zeros_like(c) for c in coeffs[1:]], WAVELET, mode="symmetric")[:len(sig)]
    rec[~valid] = np.nan
    return rec


def read_chrom_signal(bw, chrom, nb):
    arr = np.full(nb, np.nan)
    for s, e, v in bw.intervals(chrom) or []:
        i = s // BIN_BP
        if 0 <= i < nb:
            arr[i] = v
    return arr


out_path, *frac_paths = sys.argv[1:]
assert len(frac_paths) == 6, "need six fraction bigWigs in the order G1b S1 S2 S3 S4 G2"
bws = [pyBigWig.open(p) for p in frac_paths]
lens = bws[0].chroms()
chroms = [c for c in CHROMS if c in lens]
out = pyBigWig.open(out_path, "w")
out.addHeader([(c, lens[c]) for c in chroms])
for chrom in chroms:
    nb = lens[chrom] // BIN_BP
    stack = np.stack([read_chrom_signal(bw, chrom, nb) for bw in bws])
    clean = np.where(np.isfinite(stack), stack, 0.0)
    tot = clean.sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        pct = np.where(tot > 0, clean / tot * 100.0, np.nan)
    rt = np.where(tot > 0, (pct * WEIGHTS[:, None]).sum(axis=0) / 100.0, np.nan)
    sm = wavelet_smooth_J7(rt)
    ok = np.isfinite(sm)
    starts = np.arange(nb, dtype=np.int64)[ok] * BIN_BP
    if ok.any():
        out.addEntries([chrom] * int(ok.sum()), starts.tolist(), ends=(starts + BIN_BP).tolist(),
                       values=sm[ok].astype(np.float32).tolist())
out.close()
print(f"wrote {out_path}: {len(chroms)} chromosomes")
