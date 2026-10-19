# Step 9: SBS mutational signatures of the clean SNVs.
# Trinucleotide context is read from the local GATK hg38 FASTA with Rsamtools::scanFa, so the
# ~3 GB BSgenome.Hsapiens.UCSC.hg38 package is not needed. Fitting uses MutationalPatterns.
# Run from the project root:  Rscript scripts/09_signatures.R

suppressPackageStartupMessages({
  library(MutationalPatterns); library(VariantAnnotation); library(Rsamtools); library(ggplot2)
})
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# 1. SNVs only (indels have no 96-class context)
vcf <- readVcf("vcf/HCC1395.clean.vcf.gz")
vcf <- vcf[isSNV(vcf)]
gr  <- rowRanges(vcf)
cat("SNVs:", length(gr), "\n")   # 1,140

# 2. Fetch the base before and after each SNV from the reference
ctx <- as.character(scanFa(FaFile("ref/Homo_sapiens_assembly38.fasta"), resize(gr, width = 3, fix = "center")))
ref <- as.character(gr$REF)
alt <- as.character(unlist(gr$ALT))
stopifnot(all(substr(ctx, 2, 2) == ref))   # coordinates and genome build agree

# 3. Collapse to the pyrimidine (C/T) strand: G>T on + strand is the same event as C>A on - strand
comp    <- function(x) chartr("ACGT", "TGCA", x)
revcomp <- function(x) vapply(x, function(s) comp(paste(rev(strsplit(s, "")[[1]]), collapse = "")), "")
flip <- ref %in% c("G", "A")
ctx[flip] <- revcomp(ctx[flip]); ref[flip] <- comp(ref[flip]); alt[flip] <- comp(alt[flip])

# 4. Count the 96 classes (6 substitutions x 4 five-prime x 4 three-prime bases)
subs  <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
bases <- c("A", "C", "G", "T")
levels96 <- with(expand.grid(b3 = bases, b5 = bases, sub = subs), paste0(b5, "[", sub, "]", b3))
type <- paste0(substr(ctx, 1, 1), "[", ref, ">", alt, "]", substr(ctx, 3, 3))
mut_mat <- matrix(table(factor(type, levels = levels96)), ncol = 1, dimnames = list(levels96, "HCC1395"))
stopifnot(sum(mut_mat) == length(gr))
print(table(paste0(ref, ">", alt)))   # C>G and C>T dominate

ggsave("results/figures/96_profile.png", plot_96_profile(mut_mat), width = 12, height = 4, dpi = 150)
saveRDS(mut_mat, "results/tables/mut_mat.rds")

# 5. Refitting. get_known_signatures() returns a matrix without row names -> add them (same order).
signatures <- get_known_signatures()
rownames(signatures) <- levels96

# a) all 60 COSMIC SBS signatures, strict refit: over-fits (UV, platinum, haloalkane signatures appear)
fit_all <- fit_to_signatures_strict(mut_mat, signatures, max_delta = 0.004)
c_all <- fit_all$fit_res$contribution
print(round(sort(c_all[c_all[, 1] > 0, 1] / sum(c_all), decreasing = TRUE), 3))

# b) breast-relevant signatures, with vs without SBS3 (HR deficiency)
breast_sigs <- c("SBS1", "SBS2", "SBS3", "SBS5", "SBS8", "SBS13", "SBS17a", "SBS17b", "SBS18", "SBS40")
for (s in list(breast_sigs, setdiff(breast_sigs, "SBS3"))) {
  f <- fit_to_signatures(mut_mat, signatures[, s])
  cat("\n", if ("SBS3" %in% s) "with SBS3" else "without SBS3", "\n")
  print(round(sort(f$contribution[, 1] / sum(f$contribution), decreasing = TRUE), 3))
  cat("cosine similarity:", round(cos_sim(mut_mat[, 1], f$reconstructed[, 1]), 3), "\n")
}
# Result: APOBEC (SBS2 + SBS13) ~1/3 in every fit; SBS3 is replaced by flat SBS5/SBS40 with almost no
# loss of fit (0.962 vs 0.955) -> with ~1,100 SNVs, HR deficiency cannot be judged from signatures.

ggsave("results/figures/flat_vs_peaked_signatures.png",
       plot_96_profile(signatures[, c("SBS2", "SBS13", "SBS3", "SBS5", "SBS40")]),
       width = 12, height = 9, dpi = 150)
