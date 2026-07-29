suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
})

cat("======================================================================\n")
cat("PI METRIC REPORT: TOTAL CAPTURED FEATURES VS. SIGNIFICANT DEGs\n")
cat("======================================================================\n\n")

excel_file <- "Drought fractional counts.xlsx"
all_rows <- read.xlsx(excel_file, sheet = 1, colNames = FALSE)

gene_id_rows <- all_rows[, 1]
valid_gene_ids <- gene_id_rows[!is.na(gene_id_rows) & !grepl("serial|code|rep|sample|gene|id", gene_id_rows, ignore.case=TRUE)]

total_captured <- length(unique(valid_gene_ids))
cat("-> TOTAL GENES CAPTURED IN MATRIX : ", total_captured, "\n\n")

# Access downstream DEG lists from the parent directory
deg_files <- c(
  "Leaf"  = "../Supplementary_Data_Drought_LEAF.csv",
  "Root"  = "../Supplementary_Data_Drought_ROOT.csv",
  "Apex"  = "../Supplementary_Data_Drought_APEX.csv"
)

all_deg_ids <- c()
cat("-> Breakdown of Significant DEGs (padj <= 0.05) by Tissue:\n")
for (tissue in names(deg_files)) {
  file_path <- deg_files[tissue]
  if (file.exists(file_path)) {
    df <- read.csv(file_path, stringsAsFactors = FALSE)
    sig_genes <- df %>% filter(padj <= 0.05) %>% pull(Gene_ID)
    cat("   - ", tissue, " Significant DEGs: ", length(unique(sig_genes)), "\n", sep="")
    all_deg_ids <- c(all_deg_ids, sig_genes)
  }
}

if (length(all_deg_ids) > 0) {
  total_unique_degs <- length(unique(all_deg_ids))
  cat("\n-> TOTAL UNIQUE DEGs COMBINED ACROSS ALL TISSUES: ", total_unique_degs, "\n")
  cat("   (This represents ", round((total_unique_degs / total_captured) * 100, 1), 
      "% of the captured transcriptome showing a significant drought response)\n", sep="")
}
cat("======================================================================\n")
