"""Step 11: ploidy, B-allele frequency (BAF), allele-specific views and an HRD-LOH estimate.

Run from the project root:  python scripts/11_cnv_baf_hrd.py
Inputs: CNVkit output (step 10) and the GATK pileup tables at common SNPs (step 5).

Two signals are combined:
  * depth log2 ratio (tumor vs normal)  -> total copy number, relative to the tumor's average ploidy
  * BAF at SNPs heterozygous in the normal, measured inside the tumor -> balance of the two haplotypes
BAF is a within-sample ratio, so it is immune to the library-quality differences that make the depth
ratio noisy here (single normal, no panel of normals).
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIG = "results/figures/"
TAB = "results/tables/"
AUTOS = [f"chr{i}" for i in range(1, 23)]

cnr = pd.read_csv("cnv/HCC1395.dedup.cnr", sep="\t")
cns = pd.read_csv("cnv/HCC1395.dedup.cns", sep="\t")

# ---- 1. Ploidy from the spacing of segment levels -------------------------------------------
seg = cns[cns["chromosome"].isin(AUTOS)].copy()
seg["length_mb"] = (seg["end"] - seg["start"]) / 1e6
fig, ax = plt.subplots(figsize=(8, 4))
ax.hist(seg["log2"], bins=np.arange(-2, 2, 0.05), weights=seg["length_mb"], color="grey")
# Main peak at ~ -0.5 = 2 copies; then 3 copies at -0.5+log2(3/2), 4 copies at -0.5+1, 1 copy at -1.5
for n in range(1, 7):
    x = -0.5 + np.log2(n / 2)
    ax.axvline(x, color="red", ls="--", lw=1)
    ax.text(x, ax.get_ylim()[1] * 0.95, str(n), color="red", ha="center")
ax.set_xlabel("segment log2 ratio"); ax.set_ylabel("total length (Mb)")
ax.set_title("Segment levels: 2 copies at log2 ~ -0.5 -> mean ploidy ~ 2 x 2^0.5 = 2.8")
fig.tight_layout(); fig.savefig(FIG + "segment_level_histogram.png", dpi=150)

# ---- 2. BAF at germline heterozygous SNPs -----------------------------------------------------
def read_pileup(path):
    df = pd.read_csv(path, sep="\t", comment="#")
    df["depth"] = df["ref_count"] + df["alt_count"]
    df["af"] = df["alt_count"] / df["depth"]
    return df[["contig", "position", "depth", "af"]]

n = read_pileup("qc/HCC1395BL.pileups.table")
t = read_pileup("qc/HCC1395.pileups.table")
m = n.merge(t, on=["contig", "position"], suffixes=("_n", "_t"))
het = m[(m["depth_n"] >= 20) & m["af_n"].between(0.3, 0.7) & (m["depth_t"] >= 20)].copy()
print("shared sites:", len(m), " het in normal:", len(het))                      # 31,940 / 4,753
print("chrX sites:", (m["contig"] == "chrX").sum(),
      " chrX het in normal:", (het["contig"] == "chrX").sum())                     # 581 / 0 -> normal X LOH

chroms = AUTOS + ["chrX"]
h = het[het["contig"].isin(chroms)].copy()
h["ci"] = h["contig"].map({c: i for i, c in enumerate(chroms)})
h = h.sort_values(["ci", "position"])
offsets = h.groupby("ci")["position"].max().cumsum().shift(fill_value=0)
h["gpos"] = h["position"] + h["ci"].map(offsets)
fig, ax = plt.subplots(figsize=(14, 3.5))
ax.scatter(h["gpos"], h["af_t"], s=2, alpha=0.3, color="black")
for off in offsets:
    ax.axvline(off, color="grey", lw=0.5)
mids = [offsets[c] + h.loc[h["ci"] == c, "position"].max() / 2 for c in offsets.index]
ax.set_xticks(mids); ax.set_xticklabels([chroms[c].replace("chr", "") for c in offsets.index])
ax.set_ylim(0, 1); ax.set_ylabel("tumor BAF")
ax.set_title("Tumor BAF at SNPs heterozygous in the normal (points at 0/1 = LOH)")
fig.tight_layout(); fig.savefig(FIG + "baf_genome.png", dpi=150)

# ---- 3. Depth + BAF for single chromosomes ---------------------------------------------------
def zoom(chrom, marks, out):
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 5), sharex=True)
    c = cnr[cnr["chromosome"] == chrom]
    a1.scatter(c["start"] / 1e6, c["log2"], s=1, alpha=0.2, color="grey")
    for _, s in cns[cns["chromosome"] == chrom].iterrows():
        a1.plot([s["start"] / 1e6, s["end"] / 1e6], [s["log2"], s["log2"]], color="orange", lw=3)
    a1.set_ylim(-2, 2); a1.set_ylabel("log2 ratio")
    hh = het[het["contig"] == chrom]
    a2.scatter(hh["position"] / 1e6, hh["af_t"], s=4, color="black")
    a2.set_ylim(0, 1); a2.set_ylabel("tumor BAF"); a2.set_xlabel(f"{chrom} (Mb)")
    for name, pos in marks.items():
        for a in (a1, a2):
            a.axvline(pos / 1e6, color="red", ls="--")
        a1.text(pos / 1e6, 1.7, name, color="red")
    fig.tight_layout(); fig.savefig(out, dpi=150)

zoom("chr17", {"TP53": 7675088, "BRCA1": 43057078}, FIG + "chr17_depth_baf.png")
zoom("chr13", {"BRCA2": 32339132}, FIG + "chr13_depth_baf.png")

# ---- 4. HRD-LOH: LOH regions > 15 Mb that are not whole-chromosome ------------------------------
def seg_loh(r):
    hh = het[(het["contig"] == r["chromosome"]) & (het["position"] >= r["start"]) & (het["position"] < r["end"])]
    if len(hh) < 5:
        return pd.Series({"n_het": len(hh), "loh_frac": np.nan})
    return pd.Series({"n_het": len(hh), "loh_frac": ((hh["af_t"] < 0.1) | (hh["af_t"] > 0.9)).mean()})

seg[["n_het", "loh_frac"]] = seg.apply(seg_loh, axis=1)
seg["state"] = np.where(seg["loh_frac"].isna(), "unknown", np.where(seg["loh_frac"] >= 0.8, "LOH", "het"))

runs = []
for c, g in seg[seg["state"] != "unknown"].sort_values("start").groupby("chromosome"):
    cur = None
    for _, s in g.iterrows():
        if s["state"] == "LOH":
            cur = [c, s["start"], s["end"]] if cur is None else [c, cur[1], s["end"]]
        elif cur is not None:
            runs.append(cur); cur = None
    if cur is not None:
        runs.append(cur)

runs = pd.DataFrame(runs, columns=["chrom", "start", "end"])
runs["len_mb"] = (runs["end"] - runs["start"]) / 1e6
span = seg.groupby("chromosome").agg(s=("start", "min"), e=("end", "max"))
runs["chrom_frac"] = runs.apply(
    lambda r: (r["end"] - r["start"]) / (span.loc[r["chrom"], "e"] - span.loc[r["chrom"], "s"]), axis=1)
runs = runs.sort_values("len_mb", ascending=False)
runs.round(2).to_csv(TAB + "loh_regions.tsv", sep="\t", index=False)
hrd_loh = ((runs["len_mb"] > 15) & (runs["chrom_frac"] < 0.9)).sum()
print(seg["state"].value_counts())
print("HRD-LOH =", hrd_loh)   # 20 (conservative re-count ~10-14; threshold for HR deficiency ~10)
