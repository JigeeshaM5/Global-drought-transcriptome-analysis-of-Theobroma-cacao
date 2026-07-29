# 01_Metadata_Design.R
suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

cat("======================================================================\n")
cat("BASELINE: RECONSTRUCTING COHORT METADATA & PCA PATHS\n")
cat("======================================================================\n\n")

file_path <- "Drought fractional counts.xlsx"

if(!file.exists(file_path)) {
  cat("ERROR: Target file '", file_path, "' not found in the current directory.\n", sep="")
  quit(save="no", status=1)
}

# --- 1. PARSE MULTI-ROW METADATA HEADERS ---
cat("-> Extracting experimental sample attributes from multi-row sheet header...\n")
header_rows <- read.xlsx(file_path, sheet = 1, rows = 1:4, colNames = FALSE)

# Generate an accurate map of the sample attributes
sample_metadata <- data.frame(
  Column_Index = 2:ncol(header_rows),
  Code = as.character(header_rows[1, 2:ncol(header_rows)]),
  Description = as.character(header_rows[2, 2:ncol(header_rows)]),
  stringsAsFactors = FALSE
)

# Use strict regex to pull out Tissue, Time point, and Replicate identities
sample_metadata <- sample_metadata %>%
  mutate(
    Tissue = toupper(str_match(Description, "drought-([a-z]+)-")[,2]),
    TimePoint = str_match(Description, "-(T[1-7])-")[,2],
    Rep = str_match(Description, "rep #([1-3])")[,2]
  )

# Print verification check to the PI
cat("   Verification: Detected", nrow(sample_metadata), "total sample column tracks.\n")
cat("   Tissues found:", paste(unique(sample_metadata$Tissue), collapse=", "), "\n")
cat("   Time points found:", paste(unique(sample_metadata$TimePoint), collapse=", "), "\n\n")

# --- 2. LOAD NUMERIC EXPRESSION MATRIX ---
cat("-> Loading count table layers starting at row 5...\n")
counts_raw <- read.xlsx(file_path, sheet = 1, startRow = 5, colNames = FALSE)
gene_ids <- counts_raw[, 1]
counts_matrix <- as.matrix(counts_raw[, 2:ncol(counts_raw)])
colnames(counts_matrix) <- sample_metadata$Code
rownames(counts_matrix) <- gene_ids

cat("   Matrix dimensions: ", nrow(counts_matrix), " transcripts x ", ncol(counts_matrix), " samples.\n\n")

# --- 3. VARIANCE STABILIZATION & HIGH-VARIABLE GENE SELECTION ---
cat("-> Applying log2 transformation and extracting top 1,000 highly variable features...\n")
log_counts <- log2(counts_matrix + 1)
gene_variance <- apply(log_counts, 1, var)
top_features_matrix <- log_counts[order(gene_variance, decreasing = TRUE)[1:1000], ]

# --- 4. EXECUTE COHORT PRINCIPAL COMPONENT ANALYSIS ---
cat("-> Computing un-averaged cohort PCA projection coordinate values...\n")
pca_calc <- prcomp(t(top_features_matrix), scale. = TRUE)

pca_df <- data.frame(pca_calc$x[, 1:2]) %>%
  mutate(Code = rownames(.)) %>%
  inner_join(sample_metadata, by = "Code")

# Calculate centroid tracks for each time-slice to draw smooth directional paths
trajectory_paths <- pca_df %>%
  group_by(Tissue, TimePoint) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = 'drop') %>%
  arrange(Tissue, TimePoint)

# --- 5. VISUALIZE GRAPHICAL TRAJECTORY VECTOR MAP ---
cat("-> Constructing high-resolution directional PCA trajectory canvas...\n")
pc1_var <- round(summary(pca_calc)$importance[2,1] * 100, 1)
pc2_var <- round(summary(pca_calc)$importance[2,2] * 100, 1)

p_trajectory <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Tissue)) +
  # Plot individual biological replicates as subtle point clouds
  geom_point(alpha = 0.4, size = 2) +
  # Overlay continuous directional arrows representing chronological stress progression
  geom_path(data = trajectory_paths, aes(x = PC1, y = PC2, group = Tissue), linewidth = 1.3,
            arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +
  # Add chronological labels text markers to key intervals (T1 and T7)
  geom_text(data = trajectory_paths %>% filter(TimePoint %in% c("T1", "T7")),
            aes(label = TimePoint), color = "black", fontface = "bold", size = 3, vjust = -1) +
  scale_color_manual(values = c("LEAF" = "#2ecc71", "ROOT" = "#e67e22", "APEX" = "#9b59b6")) +
  labs(
    title = "Spatial-Temporal Transcriptome Trajectories in Cacao",
    subtitle = "Continuous arrows trace the chronological stress progression from baseline (T1) to severe water deficit (T7)",
    x = paste0("Principal Component 1 (", pc1_var, "% Explained Variance)"),
    y = paste0("Principal Component 2 (", pc2_var, "% Explained Variance)")
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(face = "italic", size = 9.5, color = "grey30"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

output_fig <- "01_Spatial_Temporal_PCA_Trajectories.pdf"
ggsave(output_fig, plot = p_trajectory, width = 8, height = 5.5, device = "pdf")
cat("   Success! Publication baseline saved to:", output_fig, "\n")
cat("======================================================================\n")
