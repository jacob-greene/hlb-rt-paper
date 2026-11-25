"""QC: compare the released ENCODE UW Repli-seq tracks with the tracks rebuilt by wavelet_rt.py.

Outputs in <outdir>/:
  gm12878_track.{pdf,png}       GM12878: ENCODE track (ENCFF001GOE) and the rebuilt track along chr1 and a 20 Mb window,
                                and a scatter of all 10 kb bins
  track_correlation.{pdf,png}   Pearson r of every ENCODE track against every rebuilt track (10 kb bins, chr1-22 and X)
  track_correlation.csv         the matrix
  track_correlation_pairs.csv   per track: r with its own rebuilt track, and the highest r with any other rebuilt track

Signals are read with pyBigWig (mean per 10 kb bin). The rebuilt tracks are on a 0-1 scale; they are multiplied by 100
in the plots to match the ENCODE scale. Pearson r does not depend on this scale.

Usage: qc_tracks.py <outdir> <pairs.csv>   (pairs.csv columns: label, encode, reprocessed)
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
import seaborn as sns

plt.rcParams["pdf.fonttype"] = 42
BIN = 10_000
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
outdir, pairs_csv = sys.argv[1:]
out = Path(outdir)
out.mkdir(parents=True, exist_ok=True)
pairs = pd.read_csv(pairs_csv).sort_values("label").reset_index(drop=True)


def binned(path, chroms=CHROMS):
    bw = pyBigWig.open(str(path))
    parts = []
    for c in chroms:
        n = bw.chroms()[c] // BIN
        parts.append(np.array(bw.stats(c, 0, n * BIN, type="mean", nBins=n), dtype=float))
    bw.close()
    return np.concatenate(parts)


enc = {r.label: binned(r.encode) for r in pairs.itertuples()}
rep = {r.label: binned(r.reprocessed) for r in pairs.itertuples()}
labels = list(pairs.label)
mat = pd.DataFrame(index=labels, columns=labels, dtype=float)
for a in labels:
    for b in labels:
        ok = np.isfinite(enc[a]) & np.isfinite(rep[b])
        mat.loc[a, b] = np.corrcoef(enc[a][ok], rep[b][ok])[0, 1]
mat.index.name = "encode_track"
mat.to_csv(out / "track_correlation.csv")
pd.DataFrame({"label": labels,
              "r_own": [mat.loc[l, l] for l in labels],
              "r_best_other": [mat.loc[l].drop(l).max() for l in labels],
              "best_other": [mat.loc[l].drop(l).idxmax() for l in labels]}).to_csv(out / "track_correlation_pairs.csv", index=False)

fig, ax = plt.subplots(figsize=(7.5, 6.5))
sns.heatmap(mat, cmap="viridis", vmin=0, vmax=1, square=True, ax=ax, cbar_kws={"label": "Pearson r (10 kb bins)"},
            annot=True, fmt=".2f", annot_kws={"size": 5})
ax.set_xlabel("Rebuilt from the six UCSC fraction signals (wavelet_rt.py)")
ax.set_ylabel("ENCODE released wavelet-smoothed signal")
ax.set_title("UW Repli-seq: ENCODE vs rebuilt tracks")
fig.tight_layout()
fig.savefig(out / "track_correlation.pdf")
fig.savefig(out / "track_correlation.png", dpi=300)

gm = pairs[pairs.label.str.upper() == "GM12878"].iloc[0]
e1, r1 = binned(gm.encode, ["chr1"]), binned(gm.reprocessed, ["chr1"]) * 100
x = np.arange(len(e1)) * BIN / 1e6
fig = plt.figure(figsize=(12, 6.5))
gs = fig.add_gridspec(2, 3)
for row, (lo, hi) in enumerate([(0, x.max()), (100, 120)]):
    ax = fig.add_subplot(gs[row, :2])
    m = (x >= lo) & (x <= hi)
    ax.plot(x[m], e1[m], color="#1f77b4", lw=0.8, label="ENCODE ENCFF001GOE (released)")
    ax.plot(x[m], r1[m], color="#d62728", lw=0.8, ls="--", label="rebuilt (wavelet_rt.py) x 100")
    ax.set_xlim(lo, hi)
    ax.set_ylabel("RT signal (high = early)")
    ax.set_xlabel("chr1 position (Mb)")
    ax.set_title(f"GM12878, chr1 {'whole chromosome' if row == 0 else f'{lo}-{hi} Mb'}")
    if row == 0:
        ax.legend(loc="lower left", fontsize=8)
ea, ra = enc[gm.label], rep[gm.label] * 100
ok = np.isfinite(ea) & np.isfinite(ra)
ax = fig.add_subplot(gs[:, 2])
ax.hexbin(ea[ok], ra[ok], gridsize=80, bins="log", cmap="Greys")
lim = [min(ea[ok].min(), ra[ok].min()), max(ea[ok].max(), ra[ok].max())]
ax.plot(lim, lim, color="red", lw=0.7, ls=":")
ax.set_xlabel("ENCODE released (10 kb bins)")
ax.set_ylabel("rebuilt x 100 (10 kb bins)")
ax.set_title(f"GM12878, chr1-22 and X\nPearson r = {np.corrcoef(ea[ok], ra[ok])[0, 1]:.4f}, n = {ok.sum():,} bins")
fig.tight_layout()
fig.savefig(out / "gm12878_track.pdf")
fig.savefig(out / "gm12878_track.png", dpi=300)
print(f"{len(labels)} track pairs -> {out}")
