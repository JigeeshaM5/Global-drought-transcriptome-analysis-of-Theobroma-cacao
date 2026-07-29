suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

cat("======================================================================\n")
cat("ALIGNMENT & VECTOR COERCION\n")
cat("======================================================================\n\n")

file_path <- "Drought fractional counts.xlsx"
all_rows <- read.xlsx(file_path, sheet = 1, colNames = FALSE, skipEmptyRows = FALSE)

# Detect metadata rows safely
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
  )
sample_metadata$Tissue[is.na(sample_metadata$Tissue)] <- "UNKNOWN"
sample_metadata$TimePoint[is.na(sample_metadata$TimePoint)] <- "TP"
sample_metadata$Rep[is.na(sample_metadata$Rep)] <- "1"

# Identify the exact row where numerical counts begin
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

# Coerce to strict numeric matrix
counts_numeric <- matrix(suppressWarnings(as.numeric(as.matrix(raw_mat_block))), nrow = nrow(raw_mat_block), ncol = ncol(raw_mat_block))
rownames(counts_numeric) <- gene_ids
colnames(counts_numeric) <- sample_metadata$Code

valid_rows <- apply(counts_numeric, 1, function(x) !all(is.na(x)))
counts_numeric <- counts_numeric[valid_rows, ]
counts_numeric[is.na(counts_numeric)] <- 0

log_counts <- log2(counts_numeric + 1)
gene_variance <- apply(log_counts, 1, var)
top_features_matrix <- log_counts[order(gene_variance, decreasing = TRUE)[1:1000], ]

pca_calc <- prcomp(t(top_features_matrix), scale. = TRUE)
pca_df <- data.frame(pca_calc$x[, 1:2]) %>%
  mutate(Code = rownames(.)) %>%
  inner_join(sample_metadata, by = "Code")

trajectory_paths <- pca_df %>%
  filter(Tissue != "UNKNOWN" & TimePoint != "TP") %>%
  group_by(Tissue, TimePoint) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = 'drop') %>%
  arrange(Tissue, TimePoint)

pc1_var <- round(summary(pca_calc)$importance[2,1] * 100, 1)
pc2_var <- round(summary(pca_calc)$importance[2,2] * 100, 1)

p_trajectory <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Tissue)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_path(data = trajectory_paths, aes(x = PC1, y = PC2, group = Tissue), linewidth = 1.3,
            arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +
  geom_text(data = trajectory_paths %>% filter(TimePoint %in% c("T1", "T7")),
            aes(label = TimePoint), color = "black", fontface = "bold", size = 3.5, vjust = -1.2) +
  scale_color_manual(values = c("LEAF" = "#2ecc71", "ROOT" = "#e67e22", "APEX" = "#9b59b6")) +
  labs(
    title = "Spatial-Temporal Transcriptome Trajectories in Cacao",
    subtitle = "Continuous paths map chronological drought progression across vegetative organ profiles (T1 -> T7)",
    x = paste0("Principal Component 1 (", pc1_var, "% Variance Expl.)"),
    y = paste0("Principal Component 2 (", pc2_var, "% Variance Expl.)")
  ) +
  theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold", size = 12), legend.position = "right")

ggsave("01_Spatial_Temporal_PCA_Trajectories.pdf", plot = p_trajectory, width = 8, height = 5.5, device = "pdf")
cat("-> Success! Robust PCA Figure 1 generated.\n\n")
