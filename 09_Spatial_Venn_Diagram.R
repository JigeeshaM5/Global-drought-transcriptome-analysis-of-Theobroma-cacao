# 09_Spatial_Venn_Diagram.R
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggVennDiagram)
})

cat("-> Generating High-Impact Spatial Venn Diagram (Counts Only)...\n")

leaf_genes <- if(file.exists("HighStringency_DEGs_LEAF.csv")) read.csv("HighStringency_DEGs_LEAF.csv")$Gene_ID else c()
root_genes <- if(file.exists("HighStringency_DEGs_ROOT.csv")) read.csv("HighStringency_DEGs_ROOT.csv")$Gene_ID else c()
apex_genes <- if(file.exists("HighStringency_DEGs_APEX.csv")) read.csv("HighStringency_DEGs_APEX.csv")$Gene_ID else c()

venn_list <- list(Leaf = leaf_genes, Root = root_genes, Apex = apex_genes)

# Added label = "count" to remove the misleading percentages
p_venn <- ggVennDiagram(venn_list, label_alpha = 0, edge_size = 1, set_size = 5, label = "count") +
  scale_fill_gradient(low = "#F0F0F0", high = "#009E73") +
  scale_color_manual(values = c("black", "black", "black")) +
  theme_void() +
  theme(legend.position = "none")

ggsave("03_Spatial_Venn_Diagram.pdf", plot = p_venn, width = 6, height = 6, device = "pdf")
cat("-> Success! Fixed Venn Diagram saved to: 03_Spatial_Venn_Diagram.pdf\n")
