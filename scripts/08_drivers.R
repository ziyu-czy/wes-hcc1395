# Step 8: protein-altering somatic mutations in known breast-cancer driver genes (clean set).
# Run from the project root after step 7:  Rscript scripts/08_drivers.R
# Result: TP53 p.Arg175His (VAF ~0.99) and BRCA2 p.Glu1593Ter (VAF ~0.42, wild-type copy retained).

suppressPackageStartupMessages({ library(readr); library(dplyr) })
clean <- readRDS("vep/clean_annotated.rds")

breast_drivers <- c("TP53", "PIK3CA", "GATA3", "MAP3K1", "MAP2K4", "CDH1", "PTEN", "KMT2C", "KMT2D",
                    "AKT1", "CBFB", "RUNX1", "TBX3", "NCOR1", "ERBB2", "NF1", "RB1", "ARID1A", "FOXA1",
                    "SF3B1", "BRCA1", "BRCA2", "CDKN2A", "CTCF", "ESR1", "SPEN", "PIK3R1", "KRAS", "ZFP36L1")

print(dplyr::count(clean, IMPACT))   # HIGH 50, MODERATE 314

drivers <- clean %>%
  dplyr::filter(IMPACT %in% c("HIGH", "MODERATE"), SYMBOL %in% breast_drivers) %>%
  dplyr::select(SYMBOL, Location, Consequence, HGVSp, am_class, am_pathogenicity, gnomad_af)
print(drivers, width = Inf)
write_tsv(drivers, "results/tables/somatic_drivers.tsv")

# Read support in tumor and normal (bcftools):
#   bcftools query -r chr17:7675088,chr13:32339132 -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE:%AD:%AF]\n' vcf/HCC1395.clean.vcf.gz
# Any second hit in BRCA1/BRCA2 among filtered (non-PASS) candidates? (none found)
#   bcftools query -r chr13:32315000-32400300,chr17:43044000-43126000 \
#     -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%SAMPLE:%AD:%AF]\n' vcf/HCC1395.filtered.vcf.gz
