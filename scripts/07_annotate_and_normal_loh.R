# Step 7: summarise VEP annotation of the PASS calls and detect germline leakage caused by LOH in the normal.
#
# Input: vep/HCC1395.vep.tsv.gz, from the Ensembl VEP web interface (GRCh38, GENCODE 50) run on
#        vcf/HCC1395.pass.vcf.gz with HGVS, AlphaMissense and gnomAD Exomes v4.1.1 switched on.
# Run from the project root:  Rscript scripts/07_annotate_and_normal_loh.R
#
# Story: somatic calls should almost never be common population variants. 276 PASS calls had
# gnomAD AF > 0.001, and 272 of them sat in three blocks (6p, 16q, whole chrX). There the normal showed
# 0 alt reads while the tumor showed VAF ~0.65 -> the normal cell line (HCC1395BL) lost one haplotype in
# culture, so Mutect2 mistook heterozygous germline SNPs for somatic mutations. These regions are excluded.

suppressPackageStartupMessages({ library(readr); library(dplyr) })
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Read everything as text first; readr's type guessing can silently turn sparse columns into NA.
vep <- read_tsv("vep/HCC1395.vep.tsv.gz", na = c("", "-"),
                col_types = cols(.default = col_character()))

# VEP writes one row per variant x transcript. Keep one row per variant:
# MANE Select > Ensembl canonical > any, then the most severe IMPACT.
# (Uploaded_variation is "." because Mutect2 VCFs carry no IDs, so Location + Allele is the key.)
one <- vep %>%
  dplyr::mutate(var_id = paste(Location, Allele, sep = "_"),
                priority = dplyr::case_when(!is.na(MANE_SELECT) ~ 1, CANONICAL == "YES" ~ 2, TRUE ~ 3),
                impact_rank = match(IMPACT, c("HIGH", "MODERATE", "LOW", "MODIFIER")),
                gnomad_af = suppressWarnings(as.numeric(gnomAD_exomes_AF)),   # numeric, not text
                chr = sub(":.*", "", Location),
                pos = as.integer(sub("-.*", "", sub(".*:", "", Location)))) %>%
  dplyr::group_by(var_id) %>%
  dplyr::arrange(priority, impact_rank, .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()
stopifnot(nrow(one) == 1750)

print(dplyr::count(one, IMPACT, sort = TRUE))
print(head(dplyr::count(one, Consequence, sort = TRUE), 10))

# Common population variants among "somatic" calls, and where they are
germ <- one %>% dplyr::filter(gnomad_af > 0.001)
cat("PASS calls with gnomAD AF > 0.001:", nrow(germ), "\n")
print(germ %>% dplyr::group_by(chr) %>%
        dplyr::summarise(n = dplyr::n(), from = min(pos), to = max(pos)) %>%
        dplyr::arrange(dplyr::desc(n)))
germ %>% dplyr::transmute(paste0("chr", chr), pos) %>%
  write_tsv("vep/germ_pos.tsv", col_names = FALSE)   # for: bcftools query -R vep/germ_pos.tsv ...

# Regions where the normal has LOH (from the clusters above)
regions <- tibble::tibble(chr = c("6", "16", "X"),
                          from = c(1.5e6, 46e6, 0),
                          to   = c(54e6, 90.5e6, 156.1e6))
in_reg <- one %>% dplyr::inner_join(regions, by = "chr") %>%
  dplyr::filter(pos >= from, pos <= to)
print(one %>% dplyr::mutate(in_region = var_id %in% in_reg$var_id,
                            in_gnomad = !is.na(gnomad_af) & gnomad_af > 0) %>%
        dplyr::count(in_region, in_gnomad))

# BED (0-based) of the excluded regions, used by bcftools in 07b
write_tsv(regions %>% dplyr::transmute(chrom = paste0("chr", chr), start = as.integer(from), end = as.integer(to)),
          "vcf/normal_loh.bed", col_names = FALSE)

clean <- one %>% dplyr::filter(!var_id %in% in_reg$var_id)
cat("clean set:", nrow(clean), "\n")   # 1,270
saveRDS(clean, "vep/clean_annotated.rds")
