# ----------- 1. Set Working Directory -------------
# Set the working directory to the location of the dataset
setwd("D:/NMR metabolomics")


# ------------ 2. Load Necessary R Packages -----------------
# Essential libraries for data manipulation, clustering, and visualization
library(uwot)      # UMAP dimensionality reduction
library(RANN)      # Nearest neighbors
library(readr)     # Reading data files
library(stats)     # Basic statistical functions
library(ggplot2)   # Data visualization
library(cluster)   # Clustering and silhouette analysis
library(cowplot)
library(spdep)


# ---------- 3. Load and Explore Data --------------

## ---------------------- 3.2 Load NMR Metabolomics Data ----------------------
# Read the ukb_nmr.rds file, which contains:
# - Nuclear Magnetic Resonance (NMR) metabolomics data from UK Biobank
# - 212,847 participants (rows) × 252 metabolic measurements (columns)
# - Column descriptions:
#   * f.eid: Unique participant ID in UK Biobank
#   * f.20280.0.0-f.23410.0.0: Various metabolite measurements, including:
#     - Lipoproteins
#     - Lipids
#     - Fatty acids
#     - Amino acids
#     - Other metabolic biomarkers
# - All metabolite values are numeric measurements
ukb_nmr <- readRDS("D:/NMR metabolomics/Raw_Data/ukb_nmr.rds")

# ---------------------- 4. Load and Match Disease Reporting Data ----------------------

## ---------------------- 4.1 Load Disease Reporting Data ----------------------

# Read the ukb_firstreport.rds file, which contains the following:
# - Disease first report data from UK Biobank
# - 212,847 participants (rows) × 2,261 disease reporting fields (columns)
# - Column descriptions:
#   * f.eid: Unique participant ID in UK Biobank
#   * f.130000.0.0-f.130026.0.0 etc.: Disease report information in pairs:
#     - Even numbers (f.130000.0.0): Date of first diagnosis (chr format, "YYYY-MM-DD")
#     - Odd numbers (f.130001.0.0): Source at first diagnosis.
# - Most values are NA, indicating no diagnosis for that disease
# - Each disease has two corresponding columns:
#   1. First report date (character format)
#   2. Source at first report (number format)
ukb_firstreport <- readRDS("D:/NMR metabolomics/Raw_Data/ukb_firstreport.rds")

## ---------------------- 4.2 Match Participants Between Datasets ----------------------

# Keep only participants who have both NMR and disease data
ukb_firstreport <- ukb_firstreport[ukb_firstreport$f.eid %in% ukb_nmr$f.eid, ]  # Filter disease data to keep only participants with NMR data
ukb_nmr <- ukb_nmr[ukb_nmr$f.eid %in% ukb_firstreport$f.eid, ]                  # Filter NMR data to keep only participants with disease data

## ---------------------- 4.3 Clean NMR Data ----------------------

# Create a clean NMR dataset without participant ID
# - Remove the first column (f.eid) to keep only metabolite measurements
# - Results in a matrix of just the metabolic values
nmr_data <- ukb_nmr[,-1]

# ---------------------- 5. Create a Representative Subset of Participants ----------------------

## ---------------------- 5.1 Random Sampling----------------------

# Purpose: Create a smaller, representative subset of participants for analysis
# - Reduces computational burden while maintaining population characteristics
# - Used for initial exploration and model development

# Set random seed for reproducibility
# - Ensures we get the same random sample every time we run the code
# - 2024 is used as the seed value
set.seed(2024)

# Randomly sample 50,000 participant indices
# sample(1:nrow(ukb_nmr), 50000): 
# - 1:nrow(ukb_nmr): Create sequence from 1 to total number of participants
# - Randomly select 50,000 indices from this sequence
# sort(): Order the selected indices for easier tracking
sub_set_id <- sort(sample(1:nrow(ukb_nmr), 50000))

## ---------------------- 5.2 Create and Save Subset Dataframe ----------------------

# Create a dataframe linking UK Biobank IDs with their row indices
# - ukb_id: Original UK Biobank participant IDs (from first column of ukb_nmr)
# - data_id: Row indices in our dataset (positions in the matrix)
sub_set <- data.frame(
  ukb_id = ukb_nmr[sub_set_id, 1],  # Get participant IDs for selected rows
  data_id = sub_set_id              # Store corresponding row indices
)

# Save the subset information for future reference
# - Creates CSV file mapping UK Biobank IDs to dataset row indices
# - Allows consistent subset usage across different analyses
write.csv(sub_set, "Cluster_7C_row_subset.csv")


# ---------------------- 6. Clustering Biomarkers ----------------------


## ---------------------- 6.1 Dimensionality Reduction using UMAP ----------------------

# Perform UMAP to reduce dimensionality of NMR data
umap_result <- umap(t(nmr_data),
                    n_components = 10,
                    metric = "correlation",
                    n_neighbors = 15,
                    min_dist = 0.1,
                    spread = 1,
                    seed = 24)

## ---------------------- 6.2 Calculate Distance Matrix ----------------------

# Compute distance matrix using Manhattan distance
dist_mat <- dist(umap_result, 
                 diag = TRUE,
                 upper = TRUE,
                 method = "manhattan")

## ---------------------- 6.3 Hierarchical Clustering ----------------------

# Perform hierarchical clustering
hc <- hclust(dist_mat, method = "average")

## ---------------------- 6.4 Evaluate Silhouette Scores ----------------------

# Determine optimal number of clusters (k) based on silhouette scores
si_list <- numeric(49)
for (k in 2:50) {
  clusters <- cutree(hc, k = k)
  si <- silhouette(clusters, dist_mat)
  si_list[k-1] <- mean(si[, "sil_width"])
}

## ---------------------- 6.5 Optimal Clustering and Assignments ----------------------

# Identify optimal k and perform final clustering
optimal_k <- which.max(si_list) + 1
clusters <- cutree(hc, k = optimal_k)
write.csv(clusters, "Cluster_7C_ID.csv")

## ---------------------- 6.6 Visualization (Optional) ----------------------

# Visualize clustering results if `show_plots` is set to TRUE
show_plots <- FALSE
if (show_plots) {
  par(mfrow = c(2, 2))
  
  plot(hc, main = "Hierarchical Clustering Dendrogram")
  plot(2:50, si_list, 'l', 
       main = "Silhouette Scores vs Number of Clusters",
       xlab = "Number of Clusters (k)",
       ylab = "Average Silhouette Width")
  plot(umap_result[, 1], umap_result[, 2], 
       col = clusters, cex = 1,
       main = "Clusters in UMAP Space",
       xlab = "UMAP1", ylab = "UMAP2")
  hist(clusters, 
       breaks = 1:8 - 0.5,
       main = "Distribution of Cluster Sizes",
       xlab = "Cluster", ylab = "Frequency")
  par(mfrow = c(1, 1))
}


# ---------------------- 7. Save Cluster Data ----------------------

## ---------------------- 7.1 ECDF Transformation ----------------------

# Transform metabolite data using ECDF for uniform distribution
# - apply(nmr_data, 2, rank): Calculate ranks for each metabolite (column)
# - Divide by nrow(nmr_data) to get values between 0 and 1
# - This transforms each metabolite to follow a uniform distribution
ecdf_data <- apply(nmr_data, 2, rank) / nrow(nmr_data)

## ---------------------- 7.2 Save Raw and Transformed Data ----------------------

# Process each cluster separately and save raw and ECDF-transformed data
for (i in 1:max(clusters)) {
  # Extract raw data for current cluster
  # clusters == i: Logical vector selecting metabolites in cluster i
  cluster_data <- nmr_data[, clusters == i]
  
  # Save raw data for cluster i
  csv_filename <- paste0("Clustered_raw_C_", i, ".csv")
  write.csv(cluster_data, file = csv_filename)
  
  # Extract ECDF-transformed data for current cluster
  cluster_ecdf_data <- ecdf_data[, clusters == i]
  
  # Save transformed data for cluster i
  csv_filename <- paste0("Clustered_transformed_C_", i, ".csv")
  write.csv(cluster_ecdf_data, file = csv_filename)
}

# ---------------------- ********************************** ----------------------
# ---------------------- Matlab Manifold Fitting Line ----------------------
# ---------------------- ********************************** ----------------------

# Run manifold fitting in matlab.


# ---------------------- 8. Run UMAP and Visualize for One Cluster (No color) ----------------------

# This section performs UMAP dimensionality reduction and visualization for a specific cluster.
# The visual comparison includes raw data and rank-recovered data for the cluster.

# Specify the cluster ID to process
# C_id: The ID of the cluster being analyzed and visualized
C_id <- 5
num_neighbors <- 15

## ---------------------- 8.1 Load Cluster Data ----------------------

# Load rank-transformed data for the selected cluster
# - File format: CSV without a header
rank_filename <- paste0("MF_transformed_Output_C_", C_id, ".csv")
rank_data <- read.csv(rank_filename, header = FALSE)

# Load raw metabolite data for the selected cluster
# - File format: CSV with the first column containing row indices (to be removed)
raw_filename <- paste0("Clustered_raw_C_", C_id, ".csv")
raw_data <- read.csv(raw_filename, header = TRUE)
raw_data <- raw_data[, -1]  # Remove the first column (row indices)

## ---------------------- 8.2 Recover Data from Rank Transformation ----------------------

# Reconstruct the original data distribution for the selected cluster
# - For each metabolite, the rank-transformed values are mapped back to raw quantiles
recover_data <- rank_data * 0  # Initialize an empty matrix with the same dimensions as rank_data

# Iterate through columns (metabolites) to compute quantiles
for (i in 1:ncol(rank_data)) {
  # `quantile`: Maps rank-transformed values to corresponding raw data quantiles
  recover_data[, i] <- quantile(raw_data[, i], rank_data[, i])
}

## ---------------------- 8.3 Perform UMAP for Raw and Fitted Data ----------------------

# Apply UMAP dimensionality reduction on raw and recovered data
# - `raw_data[sub_set_id, ]`: Subset the raw data using the sample subset (if applicable)
# - `n_components = 2`: Reduce to 2 dimensions for visualization
# - `metric = "correlation"`: Use correlation-based distance
# - `n_neighbors`, `min_dist`, `spread`: UMAP parameters for neighborhood size and embedding
raw_umap <- umap(raw_data[sub_set_id, ],
                 n_components = 2, metric = "correlation",
                 n_neighbors = num_neighbors, min_dist = 0.1, spread = 1, seed = 42)

rec_umap <- umap(recover_data,
                 n_components = 2, metric = "correlation",
                 n_neighbors = num_neighbors, min_dist = 0.1, spread = 1, seed = 42)

## ---------------------- 8.4 Prepare Data for Visualization ----------------------

# Convert UMAP results to data frames for easier plotting
# - `UMAP1`, `UMAP2`: Columns representing the two UMAP dimensions
raw_umap_df <- data.frame(UMAP1 = raw_umap[, 1], UMAP2 = raw_umap[, 2])
rec_umap_df <- data.frame(UMAP1 = rec_umap[, 1], UMAP2 = rec_umap[, 2])

## ---------------------- 8.5 Plot UMAP Results (Optional)  ----------------------

# Set up the plotting area for two side-by-side subplots
show_plots <- FALSE
if (show_plots) {
  par(mfrow = c(1, 2))  # 1 row, 2 columns
  

  plot(raw_umap_df$UMAP1, raw_umap_df$UMAP2,
       col = "blue", pch = 16, cex = 0.6,
       main = paste("UMAP of Raw Data (Cluster", C_id, ")"),
       xlab = "UMAP1", ylab = "UMAP2")
  grid()  # Add grid for better visual guidance
  
  # Plot UMAP for recovered data
  # - `col = "red"`: Use red color for points
  plot(rec_umap_df$UMAP1, rec_umap_df$UMAP2,
       col = "red", pch = 16, cex = 0.6,
       main = paste("UMAP of Fitted Data (Cluster", C_id, ")"),
       xlab = "UMAP1", ylab = "UMAP2")
  grid()  # Add grid for better visual guidance
  
  # Reset plotting parameters to default (single plot layout)
  par(mfrow = c(1, 1))
}


# ---------------------- 9. Run UMAP and Visualize for Multiple Clusters ----------------------
# This section performs UMAP dimensionality reduction for multiple clusters.
# The raw and rank-recovered data are processed in a loop to generate embeddings for all clusters.

## Initialize Lists to Store UMAP Results
# - `raw_umap_list`: Stores UMAP results for raw data of each cluster
# - `rec_umap_list`: Stores UMAP results for rank-recovered data of each cluster
raw_umap_list <- list()
rec_umap_list <- list()
rf_umap_list  <- list()

## Parameters for UMAP
num_neighbors <- 15
min_dist <- 0.1
spread <- 1
seed_value <- 42  # Seed for reproducibility

## Total Number of Clusters
num_clusters <- 7

## Display Progress
cat("Processing UMAP for", num_clusters, "clusters:\n")
cat("[", paste(rep(" ", num_clusters), collapse = ""), "]\r[", sep = "")

## Loop Over Clusters
# Process each cluster (e.g., clusters 1 to 7)
for (C_id in 1:num_clusters) {
  # ---------------------- 9.1 Load Cluster Data ----------------------
  
  # Load rank-transformed data for the current cluster
  rank_filename <- paste0("MF_transformed_Output_C_", C_id, ".csv")
  rank_data <- read.csv(rank_filename, header = FALSE)
  
  # Load raw metabolite data for the current cluster
  raw_filename <- paste0("Clustered_raw_C_", C_id, ".csv")
  raw_data <- read.csv(raw_filename, header = TRUE)
  raw_data <- raw_data[, -1]  # Remove the first column (row indices)
  
  # Load raw fitted data
  rf_filename <- paste0("MF_Output_C_", C_id, "_sigx10_15.csv")
  rf_data <- read.csv(rf_filename, header = FALSE)
  
  
  # ---------------------- 9.2 Recover Data from Rank Transformation ----------------------
  
  # Initialize an empty matrix for the recovered data
  recover_data <- rank_data * 0
  
  # Iterate through columns (metabolites) to reconstruct the original data distribution
  for (i in 1:ncol(rank_data)) {
    recover_data[, i] <- quantile(raw_data[, i], rank_data[, i])
  }
  
  # ---------------------- 9.3 Perform UMAP for Raw and Fitted Data ----------------------
  
  # Apply UMAP to raw data (using the specified sample subset)
  raw_umap <- umap(raw_data[sub_set_id, ],
                   n_components = 2, metric = "correlation",
                   n_neighbors = num_neighbors, min_dist = min_dist, spread = spread, seed = seed_value)
  
  # Apply UMAP to rank-recovered data
  rec_umap <- umap(recover_data,
                   n_components = 2, metric = "correlation",
                   n_neighbors = num_neighbors, min_dist = min_dist, spread = spread, seed = seed_value)
  
  # Apply UMAP to raw-fitted data
  rf_umap  <- umap(rf_data,
                   n_components = 2, metric = "correlation",
                   n_neighbors = num_neighbors, min_dist = min_dist, spread = spread, seed = seed_value)
  
  # ---------------------- 9.4 Store UMAP Results ----------------------
  
  # Save UMAP results for the current cluster into the lists
  raw_umap_list[[C_id]] <- raw_umap
  rec_umap_list[[C_id]] <- rec_umap
  rf_umap_list[[C_id]]  <- rf_umap
  
  # Update Progress
  cat("=", sep = "")
}

## Finish Progress Bar
cat("]\n")
cat("UMAP processing completed for all clusters.\n")

## Save UMAP Results
# Save both raw and recovered UMAP lists to an RDS file for later use
saveRDS(list(raw_umap_list = raw_umap_list,
             rec_umap_list = rec_umap_list,
             rf_umap_list  = rf_umap_list), 'pns_umap_lists.rds')

