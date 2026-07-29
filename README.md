# Spatiotemporal Transcriptomic Landscape of Drought Response in *Theobroma cacao*

This repository contains the complete R-based bioinformatics pipeline used to analyze the spatial and temporal transcriptomic architecture of the drought stress response in *Theobroma cacao*. 

The analytical framework isolates genuine water-deprivation-induced transcriptional signatures across three vegetative organ systems (below-ground root, orthotropic shoot apex, and foliar leaf canopy) over a seven-point diurnal time-course (T1–T7). The pipeline moves from global differential expression to temporal clustering, functional enrichment, and system-wide co-expression networks (WGCNA).

## 📊 Overview of the Analytical Pipeline

To eliminate historical threshold bias and diurnal confounding factors, the pipeline utilizes a strict multi-tier stringency threshold (FDR ≤ 0.05, |log2FC| ≥ 1.5) and matches every drought timepoint to its exact diurnal clock-matched control. 

The biological model derived from this pipeline reveals a chronologically staggered, decentralized response:
1. **Early-Stage:** Root oxidative defense and structural remodeling.
2. **Mid-Stage:** Apical meristematic quiescence via karrikin signaling.
3. **Late-Stage:** Delayed foliar senescence and photosynthetic dismantling.

---

## Repository Contents and Logical Workflow

The scripts in this repository are numbered sequentially to reflect the logical data flow of the study. Please execute them in the following order:

### Part 1: Global Transcriptomics & Differential Expression
* **`01_Data_Preprocessing_and_PCA.R`**
  * Imports raw count matrices and metadata.
  * Performs variance-stabilizing transformations (VST).
  * Generates Principal Component Analysis (PCA) plots to define global spatial topologies (Root vs. Apex vs. Leaf).
* **`02_Differential_Expression_and_Volcano.R`**
  * Executes the clock-matched differential expression analysis (DESeq2/edgeR/limma).
  * Applies the stringent FDR ≤ 0.05 and |log2FC| ≥ 1.5 thresholds.
  * Generates tissue-specific Volcano plots.
  * Calculates cross-organ intersection networks (Venn diagrams/UpSet plots).

### Part 2: Temporal Dynamics
* **`03_Temporal_Kmeans_Clustering.R`**
  * Partitions the high-stringency DEGs into distinct temporal trajectories (Clusters 1–4) across the T1–T7 continuum.
  * Plots longitudinal expression profiles.
  * **`03b_Fisher_Exact_Overlap.R`**: Calculates the statistical intersection of temporal waves across different tissues using Fisher's Exact Test to prove chronological staggering.

### Part 3: Functional Enrichment
* **`04_Functional_Enrichment_GO_ORA.R`** 
  * *Contains the corrected Over-Representation Analysis (ORA) pipeline.*
  * Maps DEGs to Gene Ontology (GO) terms to identify specialized pathways (e.g., ROS mitigation, photosynthesis shutdown).
  * Generates spatially visualized **Lollipop Plots** for the enriched GO terms.
* **`05_Functional_Enrichment_GSEA_Cnetplots.R`** 
  * *Contains the corrected Gene Set Enrichment Analysis (GSEA) pipeline.*
  * Captures coordinated sub-threshold shifts using rank-ordered transcripts.
  * Generates **Cnetplots (Gene-Concept Networks)** mapping structural nodes to highly connected signaling loci (e.g., *LHCB5*, *HSP11*).

### Part 4: Systemic Network Topologies
* **`06_WGCNA_Coexpression_Networks.R`**
  * Constructs the scale-free Weighted Gene Co-expression Network Analysis (WGCNA).
  * Identifies module eigengenes (e.g., "magenta" early root shock, "turquoise" delayed foliar response).
  * Extracts the core hub gene matrix by filtering for maximum intramodular connectivity (kME values approaching 1.0) to isolate central topological integrators (e.g., *MBF1C*, *USPAL*, *S61G1*).

---

## 💻 Prerequisites and Dependencies

The pipeline is built in R (>= 4.1.0). Ensure the following core packages are installed via CRAN or Bioconductor prior to running the scripts:

```R
# Core Bioconductor packages
BiocManager::install(c("DESeq2", "clusterProfiler", "enrichplot", "WGCNA", "ComplexHeatmap"))

# CRAN packages
install.packages(c("tidyverse", "ggplot2", "factoextra", "RColorBrewer"))
