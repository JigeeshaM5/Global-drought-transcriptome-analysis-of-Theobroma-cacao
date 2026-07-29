# 12_GSEA_Plant_Pipelines_Fixed.R
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(clusterProfiler)
  library(ggplot2)
  library(enrichplot)
  library(readxl)
})

plant_blacklist_regex <- "retina|wing|estrous|MHC|antigen presentation|p53|salivary gland|axon|neuron|muscle|sarcomere|brain|bone|lymphocyte|t-cell|b-cell|erythrocyte|testis|prostate|mammary|heart development|behavior|learning"

# ======================================================================
# PARSE & NORMALIZE NATIVE CACAO MAP FOR GSEA
# ======================================================================
excel_file <- "Drought fractional counts.xlsx"
bg_annot <- readxl::read_xlsx(excel_file, sheet = 2, skip = 2)
colnames(bg_annot) <- make.names(colnames(bg_annot))

go_map <- bg_annot %>%
  dplyr::select(Gene.ID, GO.IDs, GO.Names) %>%
  filter(!is.na(GO.IDs) & GO.IDs != "" & GO.IDs != "no GO terms") %>%
  separate_rows(GO.IDs, sep = ";\\s*") %>%
  mutate(
    Gene.ID = gsub("\\..*$", "", str_trim(Gene.ID)),
    GO.IDs = str_trim(str_remove(GO.IDs, "^[BPFC]:")), 
    GO.Names = str_trim(GO.Names)
  )

t2n_lookup <- bg_annot %>%
  dplyr::select(GO.IDs, GO.Names) %>%
  filter(!is.na(GO.IDs) & GO.IDs != "" & GO.IDs != "no GO terms") %>%
  separate_rows(c(GO.IDs, GO.Names), sep = ";\\s*") %>%
  mutate(GO.IDs = str_trim(str_remove(GO.IDs, "^[BPFC]:")), GO.Names = str_trim(str_remove(GO.Names, "^[BPFC]:"))) %>%
  distinct(GO.IDs, .keep_all = TRUE)

t2g <- go_map %>% dplyr::select(GO.IDs, Gene.ID) %>% distinct()
t2n <- t2n_lookup %>% dplyr::select(GO.IDs, GO.Names)

# ======================================================================
# GSEA RUNNER ENGINE WITH LOCUS NORMALIZATION
# ======================================================================
run_tissue_gsea <- function(data_file, tissue_name) {
  if(!file.exists(data_file)) {
    cat(sprintf("! Warning: Could not find %s. Check the file path.\n", data_file))
    return(NULL)
  }
  
  cat(sprintf("\n-> Processing Continuous GSEA for: %s\n", tissue_name))
  df <- read.csv(data_file, stringsAsFactors = FALSE)
  colnames(df) <- toupper(colnames(df))
  
  lfc_col <- grep("LOG2FOLDCHANGE|LFC", colnames(df), value = TRUE)[1]
  p_col   <- grep("PVALUE|PVAL|PADJ", colnames(df), value = TRUE)[1]
  id_col  <- grep("GENE_ID|GENE.ID|ID", colnames(df), value = TRUE)[1]
  
  # Normalize vector arrays to baseline locus names
  df <- df %>%
    mutate(
      Clean_ID = gsub("\\..*$", "", str_trim(!!sym(id_col))),
      P_CLEAN = ifelse(is.na(!!sym(p_col)) | !!sym(p_col) == 0, 1e-300, !!sym(p_col)),
      Rank_Metric = sign(!!sym(lfc_col)) * -log10(P_CLEAN)
    ) %>%
    filter(!is.na(Rank_Metric) & !is.na(Clean_ID) & Clean_ID != "") %>%
    group_by(Clean_ID) %>%
    slice_max(order_by = abs(Rank_Metric), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  gene_list <- df$Rank_Metric
  names(gene_list) <- df$Clean_ID
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  cat(sprintf("   Continuous ranked vector built with %d genes.\n", length(gene_list)))
  
  gsea_res <- GSEA(geneList = gene_list, TERM2GENE = t2g, TERM2NAME = t2n, pvalueCutoff = 0.05, pAdjustMethod = "BH", verbose = FALSE)
  
  if(is.null(gsea_res) || nrow(as.data.frame(gsea_res)) == 0) {
    cat(sprintf("   No significant GSEA pathways for %s after correction.\n", tissue_name))
    return(NULL)
  }
  
  clean_gsea_df <- as.data.frame(gsea_res) %>% filter(!grepl(plant_blacklist_regex, Description, ignore.case = TRUE))
  out_csv <- sprintf("GSEA_Results_NativeCacao_%s.csv", tissue_name)
  write.csv(clean_gsea_df, out_csv, row.names = FALSE)
  
  if(nrow(clean_gsea_df) > 0) {
    cat("   Generating publication-grade Gene-Concept Network...\n")
    gsea_res@result <- clean_gsea_df
    
    # FIXED CNETPLOT: Removed deprecated arguments
    p_cnet <- cnetplot(gsea_res, showCategory = 5, foldChange = gene_list, node_label = "category") +
      scale_color_gradient2(low = "#3182bd", mid = "#e0e0e0", high = "#de2d26", name = "Rank Metric") +
      labs(title = sprintf("Gene-Concept Network (GSEA): %s", tissue_name),
           subtitle = "Interconnected Plant Pathways Driven by Ranked Expression Trajectories") +
      theme(plot.title = element_text(face="bold", size=12))
    
    ggsave(sprintf("Figure_GSEA_Cnetplot_%s.pdf", tissue_name), plot = p_cnet, width = 11, height = 8.5, device = "pdf")
  }
  cat(sprintf("   Success! Saved CSV and Cnetplot Network for %s\n", tissue_name))
}

run_tissue_gsea("../Supplementary_Data_Drought_LEAF.csv", "LEAF")
run_tissue_gsea("../Supplementary_Data_Drought_ROOT.csv", "ROOT")
run_tissue_gsea("../Supplementary_Data_Drought_APEX.csv", "APEX")

cat("\n======================================================================\n")
cat("GLOBAL GSEA PIPELINE COMPLETE.\n")
cat("======================================================================\n")
