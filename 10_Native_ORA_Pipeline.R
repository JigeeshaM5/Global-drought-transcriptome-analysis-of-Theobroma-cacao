# 10_Native_ORA_Pipeline.R
# ======================================================================
# PARSE & NORMALIZE NATIVE CACAO BACKGROUND (PRE-FILTERED)
# ======================================================================
excel_file <- "Drought fractional counts.xlsx"
if(!file.exists(excel_file)) stop("Error: 'Drought fractional counts.xlsx' not found!")
 
cat("-> Processing Native Cacao Background from Sheet 2...\n")
bg_annot <- readxl::read_xlsx(excel_file, sheet = 2, skip = 2)
colnames(bg_annot) <- make.names(colnames(bg_annot))
 
# 1. Build the GO Name lookup and apply the lineage filter IMMEDIATELY
t2n_lookup <- bg_annot %>%
  dplyr::select(GO.IDs, GO.Names) %>%
  filter(!is.na(GO.IDs) & GO.IDs != "" & GO.IDs != "no GO terms") %>%
  separate_rows(c(GO.IDs, GO.Names), sep = ";\\s*") %>%
  mutate(
    GO.IDs = str_trim(str_remove(GO.IDs, "^[BPFC]:")), 
    GO.Names = str_trim(str_remove(GO.Names, "^[BPFC]:"))
  ) %>%
  distinct(GO.IDs, .keep_all = TRUE) %>%
  # PRE-ENRICHMENT FILTER: Purge mammalian/cross-kingdom terms
  filter(!grepl(plant_blacklist_regex, GO.Names, ignore.case = TRUE))
 
t2n <- t2n_lookup %>% dplyr::select(GO.IDs, GO.Names)

# 2. Build the GO-to-Gene map, keeping only GO IDs that survived the filter
go_map <- bg_annot %>%
  dplyr::select(Gene.ID, GO.IDs, GO.Names) %>%
  filter(!is.na(GO.IDs) & GO.IDs != "" & GO.IDs != "no GO terms") %>%
  separate_rows(GO.IDs, sep = ";\\s*") %>%
  mutate(
    Gene.ID = gsub("\\..*$", "", str_trim(Gene.ID)),
    GO.IDs = str_trim(str_remove(GO.IDs, "^[BPFC]:"))
  ) %>%
  filter(GO.IDs %in% t2n$GO.IDs) # Restrict to clean GO terms

t2g <- go_map %>% dplyr::select(GO.IDs, Gene.ID) %>% distinct()
universe_genes <- unique(t2g$Gene.ID)
 
cat(sprintf("  Background parsed: %d normalized annotations across %d unique Cacao gene loci.\n", nrow(t2g), length(universe_genes)))

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
  
  # Run enricher on the PRE-FILTERED background
  res <- enricher(gene = target_genes, universe = universe_genes, TERM2GENE = t2g, TERM2NAME = t2n, pvalueCutoff = 0.05, pAdjustMethod = "BH")
  
  if(is.null(res) || nrow(as.data.frame(res)) == 0) {
    cat(sprintf("   No terms passed statistical correction thresholds for %s.\n", tissue_name))
    return(NULL)
  }
  
  # NO FILTERING HERE. The BH FDR values are mathematically valid.
  clean_df <- as.data.frame(res) 
  out_name <- sprintf("ORA_Results_NativeCacao_%s.csv", tissue_name)
  write.csv(clean_df, out_name, row.names = FALSE)
  
  generate_lollipop_plot(
    csv_file   = out_name,
    title_text = sprintf("Native Plant GO Slim Enrichment: Top Pathways in %s", tissue_name),
    output_pdf = sprintf("Figure_Lollipop_NativeCacao_%s.pdf", tissue_name)
  )
  cat(sprintf("   Success! Saved CSV and Lollipop Figure for %s (%d plant pathways)\n", tissue_name, nrow(clean_df)))
}
