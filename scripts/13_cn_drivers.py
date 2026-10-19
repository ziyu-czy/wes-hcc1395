"""Step 13: copy number at known amplified / deleted cancer genes.

Run from the project root:  python scripts/13_cn_drivers.py
Copy number is estimated with the calibration from step 11 (log2 ~ -0.5 = 2 copies):
    CN ~ 2 * 2 ** (log2 + 0.5)
Gene coordinates are taken from the GENCODE 50 GTF (MANE Select transcript), so download
ref/gencode.v50.annotation.gtf.gz first (see README). Deep losses are then checked exon by exon
in the BAMs with scripts/13b_check_deletions.sh.
"""
import gzip
import re
import numpy as np
import pandas as pd

GENES = ["ERBB2", "MYC", "CCND1", "FGFR1", "EGFR", "PIK3CA",   # oncogenes: amplification?
         "CDKN2A", "PTEN", "RB1"]                               # tumor suppressors: deletion?

# gene spans from the MANE Select transcript records of the GTF
spans = {}
with gzip.open("ref/gencode.v50.annotation.gtf.gz", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if c[2] != "transcript" or "MANE_Select" not in c[8]:
            continue
        g = re.search(r'gene_name "([^"]+)"', c[8]).group(1)
        if g in GENES:
            spans[g] = (c[0], int(c[3]), int(c[4]))

cnr = pd.read_csv("cnv/HCC1395.dedup.cnr", sep="\t")
cns = pd.read_csv("cnv/HCC1395.dedup.cns", sep="\t")

def read_pileup(path):
    df = pd.read_csv(path, sep="\t", comment="#")
    df["depth"] = df["ref_count"] + df["alt_count"]
    df["af"] = df["alt_count"] / df["depth"]
    return df[["contig", "position", "depth", "af"]]

m = read_pileup("qc/HCC1395BL.pileups.table").merge(
    read_pileup("qc/HCC1395.pileups.table"), on=["contig", "position"], suffixes=("_n", "_t"))
het = m[(m["depth_n"] >= 20) & m["af_n"].between(0.3, 0.7) & (m["depth_t"] >= 20)]

rows = []
for g in GENES:
    c, s, e = spans[g]
    b  = cnr[(cnr["chromosome"] == c) & (cnr["end"] > s) & (cnr["start"] < e) & (cnr["gene"] != "Antitarget")]
    sg = cns[(cns["chromosome"] == c) & (cns["end"] > s) & (cns["start"] < e)]
    seg_l2 = sg["log2"].mean() if len(sg) else np.nan
    h = het[(het["contig"] == c) & (het["position"] > s - 2e6) & (het["position"] < e + 2e6)]
    loh = ((h["af_t"] < 0.1) | (h["af_t"] > 0.9)).mean() if len(h) else np.nan
    rows.append([g, c, s, e, len(b), b["log2"].median() if len(b) else np.nan,
                 seg_l2, 2 * 2 ** (seg_l2 + 0.5), len(h), loh])

res = pd.DataFrame(rows, columns=["gene", "chrom", "start", "end", "n_bins", "bin_median_log2",
                                  "segment_log2", "est_CN", "n_het_2Mb", "LOH_frac"]).round(2)
print(res.to_string(index=False))
res.to_csv("results/tables/cn_drivers.tsv", sep="\t", index=False)
# A large gap between bin_median_log2 and segment_log2 (PTEN) flags a focal change below segment resolution.
