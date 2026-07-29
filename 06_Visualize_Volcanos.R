# 06_Visualize_Volcanos.R
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

cat("======================================================================\n")
cat("GENERATING VOLCANO PLOTS (NO LABELS)\n")
cat("======================================================================\n\n")

LFC_THRESHOLD <- 1.5
FDR_THRESHOLD <- 0.05
AXIS_CAP <- 8 

process_tissue_clean <- function(file_path, tissue_name) {
  if(!file.exists(file_path)) return(NULL)
  
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  df_clean <- df %>%
    mutate(
      Tissue = tissue_name,
      LogP = -log10(padj),
      Status = case_when(
        padj <= FDR_THRESHOLD & log2FoldChange >= LFC_THRESHOLD ~ "Up-regulated",
        padj <= FDR_THRESHOLD & log2FoldChange <= -LFC_THRESHOLD ~ "Down-regulated",
        TRUE ~ "Not Significant / Low LFC"
      ),
      Plot_LFC = ifelse(log2FoldChange > AXIS_CAP, AXIS_CAP, log2FoldChange),
      Plot_LFC = ifelse(Plot_LFC < -AXIS_CAP, -AXIS_CAP, Plot_LFC)
    )
  return(df_clean)
}

cat("-> Loading data...\n")
leaf <- process_tissue_clean("../Supplementary_Data_Drought_LEAF.csv", "LEAF")
root <- process_tissue_clean("../Supplementary_Data_Drought_ROOT.csv", "ROOT")
apex <- process_tissue_clean("../Supplementary_Data_Drought_APEX.csv", "APEX")

master_df <- bind_rows(leaf, root, apex) %>%
  mutate(Tissue = factor(Tissue, levels = c("LEAF", "ROOT", "APEX")))

cat("-> Rendering clean publication canvas...\n")

p_clean <- ggplot(master_df, aes(x = Plot_LFC, y = LogP)) +
  geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD), linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_hline(yintercept = -log10(FDR_THRESHOLD), linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_point(aes(color = Status), alpha = 0.5, size = 1.0, stroke = 0) +
  facet_wrap(~ Tissue, nrow = 1) +
  scale_x_continuous(limits = c(-AXIS_CAP, AXIS_CAP), breaks = seq(-8, 8, by = 4)) +
  scale_color_manual(values = c(
    "Up-regulated" = "#D55E00",              # Okabe-Ito Vermilion
    "Down-regulated" = "#0072B2",            # Okabe-Ito Blue
    "Not Significant / Low LFC" = "grey85"   # Clean structural grey
  )) +
  labs(
    x = bquote(Log[2]~"Fold Change (Capped at ±8)"),
    y = bquote(-Log[10]~"(Adjusted P-Value)")
  ) +
  theme_classic(base_size = 14, base_family = "sans") +
  theme(
    strip.background = element_rect(fill = "#f4f4f4", color = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 13),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  ) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))

output_file <- "02_Volcanos_Uniform_Clean.pdf"
ggsave(output_file, plot = p_clean, width = 11, height = 5.0, device = "pdf")
cat("-> Success! Clean uniform plot saved to:", output_file, "\n")
