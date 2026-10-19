# Step 12b: join VEP annotation with genotypes and list variants heterozygous in the normal.
# Run from the project root:  Rscript scripts/12b_germline_hr_genes.R
# Result: BRCA1 c.5251C>T p.Arg1751Ter (ClinVar Pathogenic, expert panel), normal 86,71 -> tumor 0,67:
#         an inherited pathogenic BRCA1 allele whose wild-type partner was lost by chr17 LOH in the tumor.

suppressPackageStartupMessages({ library(readr); library(dplyr) })

vep <- read_tsv("vep/hr_genes_germline_VEP.tsv.gz", na = c("", "-"),
                col_types = cols(.default = col_character()))

one <- vep %>%
  dplyr::mutate(chrom = sub(":.*", "", Location),
                start = as.integer(sub("-.*", "", sub(".*:", "", Location))),
                priority = dplyr::case_when(!is.na(MANE_SELECT) ~ 1, CANONICAL == "YES" ~ 2, TRUE ~ 3),
                impact_rank = match(IMPACT, c("HIGH", "MODERATE", "LOW", "MODIFIER"))) %>%
  dplyr::group_by(Location, Allele) %>%
  dplyr::arrange(priority, impact_rank, .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

# VCF positions of indels point to the base before the change; VEP to the first changed base (+1)
gt <- read_tsv("vep/hr_genes_gt.tsv",
               col_names = c("chrom", "pos", "ref", "alt", "gt_n", "ad_n", "gt_t", "ad_t"),
               col_types = cols(.default = col_character())) %>%
  dplyr::mutate(chrom = sub("chr", "", chrom),
                pos = as.integer(pos),
                pos = ifelse(nchar(ref) != nchar(sub(",.*", "", alt)), pos + 1L, pos))

het_n <- one %>%
  dplyr::left_join(gt, by = c("chrom", "start" = "pos")) %>%
  dplyr::filter(gt_n == "0/1") %>%
  dplyr::select(SYMBOL, Location, Consequence, IMPACT, HGVSc, HGVSp, am_class, ClinPred, gnomAD_exomes_AF,
                ClinVar_CLNSIG, ClinVar_CLNREVSTAT, dplyr::starts_with("SpliceAI_pred_DS"),
                gt_n, ad_n, gt_t, ad_t)
print(het_n, n = Inf, width = Inf)
write_tsv(het_n, "results/tables/germline_hr_het_variants.tsv")
