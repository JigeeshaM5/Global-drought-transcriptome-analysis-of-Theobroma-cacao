# 08_Spatial_Table.R
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
})

cat("======================================================================\n")
cat("GENERATING MASTER SPATIAL DEG ANNOTATION TABLE\n")
cat("======================================================================\n\n")

load_deg <- function(file_path, prefix) {
  if(!file.exists(file_path)) return(data.frame(Gene_ID = character()))
  
  df <- read.csv(file_path, stringsAsFactors = FALSE) %>%
    select(Gene_ID, Clean_Symbol, Description, log2FoldChange, padj, Status)
  
  # Rename columns to reflect the tissue
  colnames(df)[4:6] <- paste0(prefix, c("_LFC", "_padj", "_Status"))
  return(df)
}

# Load the stringent subset files we generated in the last step
leaf <- load_deg("HighStringency_DEGs_LEAF.csv", "LEAF")
root <- load_deg("HighStringency_DEGs_ROOT.csv", "ROOT")
apex <- load_deg("HighStringency_DEGs_APEX.csv", "APEX")

# Full outer join to capture every unique DEG across the plant
master_df <- leaf %>%
  full_join(root, by = c("Gene_ID", "Clean_Symbol", "Description")) %>%
  full_join(apex, by = c("Gene_ID", "Clean_Symbol", "Description"))

# Categorize the spatial compartments logically
master_df <- master_df %>%
  mutate(
    In_Leaf = !is.na(LEAF_Status),
    In_Root = !is.na(ROOT_Status),
    In_Apex = !is.na(APEX_Status),
    
    Spatial_Compartment = case_when(
      In_Leaf & In_Root & In_Apex ~ "Universal Core",
      In_Leaf & In_Root & !In_Apex ~ "Leaf + Root Shared",
      In_Leaf & !In_Root & In_Apex ~ "Leaf + Apex Shared",
      !In_Leaf & In_Root & In_Apex ~ "Root + Apex Shared",
      In_Leaf & !In_Root & !In_Apex ~ "Leaf Unique",
      !In_Leaf & In_Root & !In_Apex ~ "Root Unique",
      !In_Leaf & !In_Root & In_Apex ~ "Apex Unique",
      TRUE ~ "Unknown"
    )
  ) %>%
  select(Gene_ID, Clean_Symbol, Spatial_Compartment, Description, 
         LEAF_LFC, LEAF_padj, LEAF_Status, 
         ROOT_LFC, ROOT_padj, ROOT_Status, 
         APEX_LFC, APEX_padj, APEX_Status) %>%
  arrange(Spatial_Compartment, Gene_ID)

# Save the final pristine table
output_csv <- "Table_S1_Master_Spatial_Stringent_DEGs.csv"
write.csv(master_df, output_csv, row.names = FALSE)
cat("-> Success! Master table saved to:", output_csv, "\n")
cat("   Total Stringent Unique DEGs captured:", nrow(master_df), "\n")
cat("======================================================================\n")
