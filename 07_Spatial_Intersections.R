# 07_Spatial_Intersections.R
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
})

cat("======================================================================\n")
cat("GENERATING INTEGRATED ACCESSIBLE INTERSECTION ANALYSIS (VENN & UPSET)\n")
cat("======================================================================\n\n")

LFC_THRESHOLD <- 1.5
FDR_THRESHOLD <- 0.05

load_stringent_genes <- function(file_path) {
  if(!file.exists(file_path)) return(character(0))
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  degs <- df %>% 
    filter(padj <= FDR_THRESHOLD & abs(log2FoldChange) >= LFC_THRESHOLD) %>%
    pull(Gene_ID)
  return(unique(degs))
}

# Load pure gene sets
leaf_genes <- load_stringent_genes("../Supplementary_Data_Drought_LEAF.csv")
root_genes <- load_stringent_genes("../Supplementary_Data_Drought_ROOT.csv")
apex_genes <- load_stringent_genes("../Supplementary_Data_Drought_APEX.csv")

# Print ground-truth calculations to terminal
cat(sprintf("Stringent DEGs captured:\n  Leaf: %d\n  Root: %d\n  Apex: %d\n\n", 
            length(leaf_genes), length(root_genes), length(apex_genes)))

# Solve explicit intersection matrices
all_universe <- unique(c(leaf_genes, root_genes, apex_genes))

# Structural categories
L_only <- length(setdiff(leaf_genes, c(root_genes, apex_genes)))
R_only <- length(setdiff(root_genes, c(leaf_genes, apex_genes)))
A_only <- length(setdiff(apex_genes, c(leaf_genes, root_genes)))

LR_shared <- length(setdiff(intersect(leaf_genes, root_genes), apex_genes))
LA_shared <- length(setdiff(intersect(leaf_genes, apex_genes), root_genes))
RA_shared <- length(setdiff(intersect(root_genes, apex_genes), leaf_genes))

Universal_Core <- length(intersect(intersect(leaf_genes, root_genes), apex_genes))

cat("---------------------------------------------------\n")
cat("INTERSECTION MATH VERIFICATION:\n")
cat("  Leaf Unique Specific:     ", L_only, "\n")
cat("  Root Unique Specific:     ", R_only, "\n")
cat("  Apex Unique Specific:     ", A_only, "\n")
cat("  Leaf-Root Shared Excl:    ", LR_shared, "\n")
cat("  Leaf-Apex Shared Excl:    ", LA_shared, "\n")
cat("  Root-Apex Shared Excl:    ", RA_shared, "\n")
cat("  UNIVERSAL CORE SHARED (3):", Universal_Core, "\n")
cat("---------------------------------------------------\n\n")

# 1. BUILD HIGH-IMPACT UPSET DATA MATRIX
upset_df <- data.frame(
  Size = c(L_only, R_only, A_only, LA_shared, LR_shared, RA_shared, Universal_Core),
  Combination = c("Leaf Only", "Root Only", "Apex Only", "Leaf+Apex", "Leaf+Root", "Root+Apex", "Universal Core")
) %>% arrange(desc(Size))

upset_df$Combination <- factor(upset_df$Combination, levels = upset_df$Combination)

p_upset <- ggplot(upset_df, aes(x = Combination, y = Size)) +
  geom_col(fill = "#009E73", width = 0.6) + # High visibility solid bar fill
  geom_text(aes(label = scales::comma(Size)), vjust = -0.5, size = 4.5, fontface = "bold") +
  labs(x = "Organ Convergence Compartments", y = "Differentially Expressed Gene Count") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5)
  )

# 2. SAVE SEPARATE CSV CORES FOR DOWNSTREAM ANALYSIS
write.csv(data.frame(Gene_ID = intersect(intersect(leaf_genes, root_genes), apex_genes)), "Universal_Core_Shared_DEGs.csv", row.names=FALSE)
write.csv(data.frame(Gene_ID = setdiff(leaf_genes, c(root_genes, apex_genes))), "Leaf_Specific_DEGs.csv", row.names=FALSE)
write.csv(data.frame(Gene_ID = setdiff(root_genes, c(leaf_genes, apex_genes))), "Root_Specific_DEGs.csv", row.names=FALSE)
write.csv(data.frame(Gene_ID = setdiff(apex_genes, c(leaf_genes, root_genes))), "Apex_Specific_DEGs.csv", row.names=FALSE)

ggsave("03_Spatial_DEG_Intersections_UpSet.pdf", plot = p_upset, width = 8, height = 5.5, device = "pdf")
cat("-> Success! High-impact intersection plot saved to: 03_Spatial_DEG_Intersections_UpSet.pdf\n")
cat("======================================================================\n")
