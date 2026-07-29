# 05_PCA_Analysis.R
suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

cat("-> Generating High-Impact, Accessible PCA Visualization...\n")

file_path <- "Drought fractional counts.xlsx"
all_rows <- read.xlsx(file_path, sheet = 1, colNames = FALSE, skipEmptyRows = FALSE)

# 1. Parse Metadata safely
sample_code_idx <- which(apply(all_rows, 1, function(r) any(grepl("Sample Code Name", r, ignore.case=TRUE))))
descriptive_idx <- which(apply(all_rows, 1, function(r) any(grepl("Sample Descriptive Name", r, ignore.case=TRUE))))
if(length(sample_code_idx) == 0) sample_code_idx <- 1
if(length(descriptive_idx) == 0) descriptive_idx <- 2

codes <- as.character(all_rows[sample_code_idx, -1])
descriptions <- as.character(all_rows[descriptive_idx, -1])

sample_metadata <- data.frame(Code = codes, Description = descriptions, stringsAsFactors = FALSE) %>%
  mutate(
    Tissue = toupper(str_match(Description, "drought-([a-zA-Z]+)-")[,2]),
    TimePoint = str_match(Description, "-(T[1-7])-")[,2],
    Rep = str_match(Description, "rep #([1-3])")[,2]
  ) %>%
  mutate(
    Time_Numeric = as.numeric(factor(TimePoint, levels = paste0("T", 1:7)))
  )

# 2. Extract and strict-coerce Numeric Matrix
numeric_start_idx <- 1
for(i in 1:nrow(all_rows)) {
  val <- all_rows[i, 2]
  if(!is.na(val) && !is.na(suppressWarnings(as.numeric(val))) && !grepl("serial|code|rep|sample", val, ignore.case=TRUE)) {
    numeric_start_idx <- i
    break
  }
}

data_rows <- all_rows[numeric_start_idx:nrow(all_rows), ]
gene_ids <- as.character(data_rows[, 1])
raw_mat_block <- data_rows[, 2:ncol(data_rows)]

counts_numeric <- matrix(suppressWarnings(as.numeric(as.matrix(raw_mat_block))), nrow = nrow(raw_mat_block), ncol = ncol(raw_mat_block))
rownames(counts_numeric) <- gene_ids
colnames(counts_numeric) <- sample_metadata$Code

valid_rows <- apply(counts_numeric, 1, function(x) !all(is.na(x)))
counts_numeric <- counts_numeric[valid_rows, ]
counts_numeric[is.na(counts_numeric)] <- 0

# 3. Log-transform and run PCA
log_counts <- log2(counts_numeric + 1)
gene_variance <- apply(log_counts, 1, var)
top_features_matrix <- log_counts[order(gene_variance, decreasing = TRUE)[1:1000], ]

pca_calc <- prcomp(t(top_features_matrix), scale. = TRUE)
pca_df <- data.frame(pca_calc$x[, 1:2]) %>%
  mutate(Code = rownames(.)) %>%
  inner_join(sample_metadata, by = "Code") %>%
  filter(!is.na(Tissue) & Tissue != "NA")

# Calculate centroid pathways
trajectory_paths <- pca_df %>%
  group_by(Tissue, Time_Numeric) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = 'drop') %>%
  arrange(Tissue, Time_Numeric)

# Variance explained for axis labels
pc1_var <- round(summary(pca_calc)$importance[2,1] * 100, 1)
pc2_var <- round(summary(pca_calc)$importance[2,2] * 100, 1)

# 4. Plot using Okabe-Ito Palette, no text in plot, mapped opacity
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  # Draw sleek centroid connecting paths (no arrows)
  geom_path(data = trajectory_paths, aes(color = Tissue), linewidth = 1.2, alpha = 0.5, show.legend = FALSE) +
  # Points with white outlines for high contrast, fill tied to tissue, opacity tied to time
  geom_point(aes(fill = Tissue, alpha = Time_Numeric), shape = 21, size = 4, color = "white", stroke = 0.6) +
  # Okabe-Ito Colorblind-Safe Palette
  scale_fill_manual(values = c("LEAF" = "#009E73", "ROOT" = "#D55E00", "APEX" = "#0072B2")) +
  scale_color_manual(values = c("LEAF" = "#009E73", "ROOT" = "#D55E00", "APEX" = "#0072B2")) +
  scale_alpha_continuous(range = c(0.2, 1), breaks = 1:7, labels = paste0("T", 1:7), name = "Time Point") +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  ) +
  # High-impact classic theme (minimalist)
  theme_classic(base_size = 14, base_family = "sans") +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold", margin = margin(t = 10, r = 10)),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black")
  ) +
  # Customize legend guides for a cleaner look
  guides(
    fill = guide_legend(override.aes = list(alpha = 1, size = 4)),
    alpha = guide_legend(override.aes = list(fill = "grey30"))
  )

ggsave("01_Spatial_Temporal_PCA_HighImpact.pdf", plot = p_pca, width = 7.5, height = 5.5, device = "pdf")
cat("-> Success! High-Impact PCA saved to: 01_Spatial_Temporal_PCA_HighImpact.pdf\n")
