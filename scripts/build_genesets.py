"""Build the protein-coding gene table and the replication-dependent histone set.

Input: GENCODE v44 basic annotation lifted to GRCh37 (gencode.v44lift37.basic.annotation.gtf.gz).
Output (in OUTDIR): all_coding.parquet (one row per protein-coding gene; canonical TSS),
histones.parquet (replication-dependent histone genes in the four clusters HIST1-HIST4).
Coordinates are kept exactly as in the GTF (1-based start), as in the original analysis.
Usage: build_genesets.py <gtf.gz> <outdir>
"""
import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

GTF = Path(sys.argv[1])
OUT = Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)

ATTR_PAT = re.compile(r'(\w+) "([^"]*)"')


def parse_attr(s):
    return dict(ATTR_PAT.findall(s))


def read_gtf(path):
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, src, ftype, start, end, score, strand, frame, attrs = parts
            if ftype not in ("transcript", "gene"):
                continue
            attr = parse_attr(attrs)
            rows.append((seqid, ftype, int(start), int(end), strand,
                         attr.get("gene_id"), attr.get("gene_name"),
                         attr.get("gene_type"), attr.get("transcript_id"),
                         attr.get("transcript_type"),
                         attr.get("transcript_support_level"),
                         attr.get("tag", "")))
    return pd.DataFrame(rows, columns=[
        "chrom", "ftype", "start", "end", "strand", "gene_id", "gene_name",
        "gene_type", "transcript_id", "transcript_type", "tsl", "tag"])


def main():
    df = read_gtf(GTF)
    print(f"  raw rows: {len(df)}")

    genes = df[df["ftype"] == "gene"].copy()
    tx = df[df["ftype"] == "transcript"].copy()
    print(f"  gene rows: {len(genes)}; transcript rows: {len(tx)}")

    # Protein-coding only
    pc_genes = genes[genes["gene_type"] == "protein_coding"].copy()
    pc_tx = tx[tx["gene_type"] == "protein_coding"].copy()
    print(f"  protein-coding gene rows: {len(pc_genes)}; tx rows: {len(pc_tx)}")

    # Pick canonical transcript per gene: prefer one tagged Ensembl_canonical
    # if present, then MANE_Select, else longest tx.
    pc_tx["length"] = pc_tx["end"] - pc_tx["start"]

    def canonical_tx(group):
        if (group["tag"].str.contains("Ensembl_canonical", na=False)).any():
            sub = group[group["tag"].str.contains("Ensembl_canonical", na=False)]
        elif (group["tag"].str.contains("MANE_Select", na=False)).any():
            sub = group[group["tag"].str.contains("MANE_Select", na=False)]
        else:
            sub = group
        return sub.sort_values("length", ascending=False).iloc[0]

    print("  selecting canonical transcript per gene…")
    # Pick canonical row indices manually to avoid groupby+apply quirks
    pc_tx["_has_canonical"] = pc_tx["tag"].str.contains("Ensembl_canonical", na=False)
    pc_tx["_has_mane"] = pc_tx["tag"].str.contains("MANE_Select", na=False)
    pc_tx["_priority"] = (
        pc_tx["_has_canonical"].astype(int) * 100 +
        pc_tx["_has_mane"].astype(int) * 10
    )
    pc_tx = pc_tx.sort_values(
        ["gene_id", "_priority", "length"], ascending=[True, False, False])
    canonical = pc_tx.drop_duplicates("gene_id", keep="first").reset_index(drop=True)
    print(f"  canonical tx rows: {len(canonical)}")

    # Build per-gene table with canonical TSS, gene body span (gene record),
    # and strand.
    canonical_short = canonical[["gene_id", "transcript_id", "tag", "tsl",
                                 "length"]].rename(columns={"length": "tx_length",
                                                            "tag": "tx_tag",
                                                            "tsl": "tx_tsl"})
    pc = pc_genes[["chrom", "start", "end", "strand", "gene_id",
                   "gene_name"]].merge(canonical_short, on="gene_id",
                                       how="left")
    pc["tss"] = np.where(pc["strand"] == "+", pc["start"], pc["end"])
    pc["gene_length"] = pc["end"] - pc["start"]
    # Strip versioning from gene_id for downstream merges
    pc["gene_id_unversioned"] = pc["gene_id"].str.split(".").str[0]
    print(f"  protein-coding gene-level rows: {len(pc)}")

    # Drop duplicates (lift37 occasionally duplicates a gene on chr_alt)
    pc = pc.sort_values(["chrom", "start"]).drop_duplicates("gene_id", keep="first")
    pc.to_parquet(OUT / "all_coding.parquet", index=False)
    pc.to_csv(OUT / "all_coding.csv", index=False)
    print(f"  saved all_coding ({len(pc)} genes)")

    # ---- Replication-dependent histone gene set ----
    # Four canonical clusters per HGNC histone-genes nomenclature (2019):
    #   HIST1  chr6p21.3-22.2 (~chr6:25.7–27.9 Mb on hg19)
    #   HIST2  chr1q21.2     (~chr1:149.7–149.86 Mb)
    #   HIST3  chr1q42.13    (~chr1:228.6–228.65 Mb)
    #   HIST4  chr12p12.3    (~chr12:14.92 Mb, single gene H4C16)
    # Symbols (post-2019 HGNC): H1-{1..6}, H2AC{1..21,25}, H2BC{1..26},
    # H3C{1..15}, H4C{1..16}, H3-4 (HIST3H3); legacy: HIST1H4B etc.
    # We INCLUDE replication-dependent histones; we EXCLUDE replication-
    # independent variants (H2AFX/H2AX, H3F3A/H3F3B = H3-3A/B, MACROH2A1/2,
    # H2AZ1/2, CENPA, H3-3, H3-5, H2AB1/2/3 (variant), H2AJ, H1-0, H1-7,
    # H1-8, H1-10 — somatic-tissue / replacement / testis variants).
    rd_pattern = re.compile(
        r"^(HIST[1-4]H[1234][A-Z]?\d*[A-Z]?|"            # legacy HIST1H4A etc.
        r"H1[-_]?[1-6]|"                                    # H1-1..H1-6 (rep-dep)
        r"H2A[CB]\d+|H2BC\d+|H2BU\d+|"                      # H2AC, H2BC, H2BU
        r"H3C\d+|H4C\d+|"                                    # H3C, H4C
        r"H3-4)$"                                             # HIST3H3 (rep-dep)
    )
    # Variant / replacement / testis-specific histones excluded even if
    # symbol pattern accidentally matches.
    excluded_variants = {"H2AX", "H2AFX", "H3F3A", "H3F3B", "H3-3A", "H3-3B",
                          "H3-5", "H3-7",  # H3 variants (rep-independent)
                          "MACROH2A1", "MACROH2A2", "H2AFZ", "H2AZ1", "H2AZ2",
                          "H2A.X", "H2A.Z", "CENPA",
                          "H2AJ", "H2AB1", "H2AB2", "H2AB3",  # variant H2A
                          "H1-0", "H1-7", "H1-8", "H1-10"}     # H1 variants

    pc["is_histone_symbol"] = (
        pc["gene_name"].str.match(rd_pattern, na=False) &
        ~pc["gene_name"].isin(excluded_variants)
    )
    # Spatial filter + per-cluster label
    in_hist1 = ((pc["chrom"] == "chr6") &
                (pc["start"] >= 25_000_000) &
                (pc["end"] <= 28_500_000))
    in_hist2 = ((pc["chrom"] == "chr1") &
                (pc["start"] >= 149_300_000) &
                (pc["end"] <= 150_000_000))
    in_hist3 = ((pc["chrom"] == "chr1") &
                (pc["start"] >= 228_400_000) &
                (pc["end"] <= 228_800_000))
    in_hist4 = ((pc["chrom"] == "chr12") &
                (pc["start"] >= 14_800_000) &
                (pc["end"] <= 15_000_000))
    pc["in_histone_cluster"] = in_hist1 | in_hist2 | in_hist3 | in_hist4
    pc["is_histone"] = pc["is_histone_symbol"] & pc["in_histone_cluster"]
    pc["histone_cluster"] = np.where(in_hist1, "HIST1",
                              np.where(in_hist2, "HIST2",
                                np.where(in_hist3, "HIST3",
                                  np.where(in_hist4, "HIST4", ""))))
    pc["histone_cluster"] = np.where(pc["is_histone"], pc["histone_cluster"], "")

    histones = pc[pc["is_histone"]].copy()
    histones.to_parquet(OUT / "histones.parquet", index=False)
    histones.to_csv(OUT / "histones.csv", index=False)
    print(f"  histones (rep-dep): {len(histones)}")
    if len(histones) < 30:
        print("  WARNING: fewer histones than expected (≥40)")
    print("  histone breakdown by cluster:")
    print("   ", histones["histone_cluster"].value_counts().to_dict())
    for cl in ["HIST1", "HIST2", "HIST3", "HIST4"]:
        sub = histones[histones["histone_cluster"] == cl]
        print(f"  {cl}: {sorted(sub['gene_name'].tolist())}")

if __name__ == "__main__":
    main()
