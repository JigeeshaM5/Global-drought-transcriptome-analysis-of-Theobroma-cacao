# 11_GSEA_Pipeline.R
# ======================================================================
# PARSE & NORMALIZE NATIVE CACAO MAP FOR GSEA (PRE-FILTERED)
# ======================================================================
excel_file <- "Drought fractional counts.xlsx"
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
  # PRE-ENRICHMENT FILTER: Purge mammalian terms from the GO Universe
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

# ======================================================================
# GSEA RUNNER ENGINE WITH LOCUS NORMALIZATION
# ======================================================================
# (Inside your run_tissue_gsea function, remove the post-filter)

  gsea_res <- GSEA(geneList = gene_list, TERM2GENE = t2g, TERM2NAME = t2n, pvalueCutoff = 0.05, pAdjustMethod = "BH", verbose = FALSE)
  
  if(is.null(gsea_res) || nrow(as.data.frame(gsea_res)) == 0) {
    cat(sprintf("   No significant GSEA pathways for %s after correction.\n", tissue_name))
    return(NULL)
  }
  
  # No filtering needed here anymore. The BH FDR values are now mathematically valid.
  clean_gsea_df <- as.data.frame(gsea_res) 
  out_csv <- sprintf("GSEA_Results_NativeCacao_%s.csv", tissue_name)
  write.csv(clean_gsea_df, out_csv, row.names = FALSE)
