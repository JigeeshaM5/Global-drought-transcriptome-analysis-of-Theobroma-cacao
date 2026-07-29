# 11_Native_ORA_Slim_Pipeline.R
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "http://cran.us.r-project.org", quiet = TRUE)
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler", update = FALSE, quiet = TRUE)
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl", repos = "http://cran.us.r-project.org", quiet = TRUE)
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos = "http://cran.us.r-project.org", quiet = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(clusterProfiler)
  library(readxl)
  library(ggplot2)
})

# Cross-kingdom contaminant purge
plant_blacklist_regex <- "retina|wing|estrous|MHC|antigen presentation|p53|salivary gland|axon|neuron|muscle|sarcomere|brain|bone|lymphocyte|t-cell|b-cell|erythrocyte|testis|prostate|mammary|heart development|behavior|learning"

# ======================================================================
# LOLLIPOP PLOT VISUALIZATION ENGINE
# ======================================================================
generate_lollipop_plot <- function(csv_file, title_text, output_pdf) {
  if(!file.exists(csv_file)) return(NULL)
  df <- read.csv(csv_file, stringsAsFactors = FALSE)
  if(nrow(df) == 0) return(NULL)
  
  plot_df <- df %>%
    arrange(p.adjust) %>%
    slice_head(n = 15) %>%
    mutate(
      Ratio_Num = sapply(GeneRatio, function(x) {
        nums <- as.numeric(strsplit(x, "/")[[1]])
        nums[1] / nums[2]
      }),
      Wrapped_Desc = stringr::str_wrap(Description, width = 45)
    ) %>%
    arrange(Ratio_Num) %>%
    mutate(Wrapped_Desc = factor(Wrapped_Desc, levels = unique(Wrapped_Desc)))
  
  p <- ggplot(plot_df, aes(x = Ratio_Num, y = Wrapped_Desc)) +
    geom_segment(aes(x = 0, xend = Ratio_Num, y = Wrapped_Desc, yend = Wrapped_Desc), 
                 color = "grey80", linewidth = 0.75) +
    geom_point(aes(size = Count, color = p.adjust), alpha = 0.95) +
    scale_color_gradientn(colors = c("#de2d26", "#fc9272", "#fee0d2"), 
                          name = "Adj. p-value",
                          limits = c(0, max(plot_df$p.adjust, 0.05))) +
    scale_size_continuous(range = c(4, 9), name = "Gene Count") +
    theme_minimal(base_size = 11) +
    labs(title = stringr::str_wrap(title_text, width = 55),
         x = "Gene Ratio (Enriched DEGs / Total Input DEGs)", y = NULL) +
    theme(
      plot.title = element_text(face = "bold", size = 11, color = "black"),
      axis.text.y = element_text(size = 9.5, color = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "right"
    )
  
  ggsave(output_pdf, plot = p, width = 8.5, height = 6.5, device = "pdf")
}

# ======================================================================
# PARSE & NORMALIZE NATIVE CACAO BACKGROUND
# ======================================================================
excel_file <- "Drought fractional counts.xlsx"
if(!file.exists(excel_file)) stop("Error: 'Drought fractional counts.xlsx' not found!")

cat("-> Processing Native Cacao Background from Sheet 2...\n")
bg_annot <- readxl::read_xlsx(excel_file, sheet = 2, skip = 2)
colnames(bg_annot) <- make.names(colnames(bg_annot))

# Normalize background IDs by stripping transcript suffixes (e.g., .1)
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
universe_genes <- unique(go_map$Gene.ID)

cat(sprintf("   Background parsed: %d normalized annotations across %d unique Cacao gene loci.\n", nrow(t2g), length(universe_genes)))

# ======================================================================
# CORE ORA RUNNER WITH ISOFORM NORMALIZATION
# ======================================================================
run_native_ora <- function(data_file, tissue_name) {
  if(!file.exists(data_file)) return(NULL)
  df <- read.csv(data_file, stringsAsFactors = FALSE)
  colnames(df) <- make.names(colnames(df))
  id_col <- grep("Gene.ID|Gene_ID|ID", colnames(df), value = TRUE)[1]
  
  # Strip transcript isoform suffixes from foreground targets to match background loci
  raw_target_genes <- df[[id_col]]
  target_genes <- unique(gsub("\\..*$", "", str_trim(raw_target_genes)))
  
  cat(sprintf("-> Running Native ORA for %s...\n", tissue_name))
  
  # Diagnostic Match Verification
  overlap_count <- sum(target_genes %in% universe_genes)
  cat(sprintf("   [Diagnostic] Locus Normalization Success: %d out of %d DEGs successfully mapped to background.\n", 
              overlap_count, length(target_genes)))
  
  if(overlap_count == 0) {
    cat("   Stopping execution: Zero mapping overlap detected. Verify file schemas.\n")
    return(NULL)
  }
  
  res <- enricher(gene = target_genes, universe = universe_genes, TERM2GENE = t2g, TERM2NAME = t2n, pvalueCutoff = 0.05, pAdjustMethod = "BH")
  if(is.null(res) || nrow(as.data.frame(res)) == 0) {
    cat(sprintf("   No terms passed statistical correction thresholds for %s.\n", tissue_name))
    return(NULL)
  }
  
  clean_df <- as.data.frame(res) %>% filter(!grepl(plant_blacklist_regex, Description, ignore.case = TRUE))
  out_name <- sprintf("ORA_Results_NativeCacao_%s.csv", tissue_name)
  write.csv(clean_df, out_name, row.names = FALSE)
  
  generate_lollipop_plot(
    csv_file   = out_name,
    title_text = sprintf("Native Plant GO Slim Enrichment: Top Pathways in %s", tissue_name),
    output_pdf = sprintf("Figure_Lollipop_NativeCacao_%s.pdf", tissue_name)
  )
  cat(sprintf("   Success! Saved CSV and Lollipop Figure for %s (%d plant pathways)\n", tissue_name, nrow(clean_df)))
}

run_native_ora("HighStringency_DEGs_LEAF.csv", "LEAF")
run_native_ora("HighStringency_DEGs_ROOT.csv", "ROOT")
run_native_ora("HighStringency_DEGs_APEX.csv", "APEX")
cat("\nNATIVE ORA PIPELINE COMPLETE.\n")
