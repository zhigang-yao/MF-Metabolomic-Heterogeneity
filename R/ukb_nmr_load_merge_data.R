setwd("D:/NMR metabolomics")

#---------------------------
# Load packages
#---------------------------

# Data Processing and Analysis
library(uwot)        # UMAP dimensionality reduction
library(RANN)        # Fast nearest neighbor search
library(readr)       # Reading and writing data files efficiently
library(stats)       # Basic statistical functions and utilities
library(dplyr)       # Data manipulation and transformation
library(tidyr)       # Data tidying and reshaping

# Clustering
library(cluster)     # Various clustering algorithms and silhouette analysis
library(dbscan)      # Density-based clustering

# Visualization
library(ggplot2)     # Grammar of graphics plotting
library(ggrepel)     # Text and label repulsion for ggplot2
library(cowplot)     # Publication-quality figure arrangement
library(gridExtra)   # Arranging multiple grid-based plots
library(scales)      # Scale functions for visualization

# Survival Analysis
library(survival)    # Survival analysis
library(forestplot)  # Create forest plots for meta-analyses


#---------------------------
# Load data
#---------------------------

ukb_all <- readRDS("D:/NMR metabolomics/Raw_Data/UKB_All_Data.rds")

Participant_ID_Subset_50K <- read.csv("D:/NMR metabolomics/Raw_Data/Cluster_7C_row_subset.csv")

ICD10_FieldID <- read.csv("D:/NMR metabolomics/Raw_Data/ICD10_FieldID.csv")

ukb_firstreport <- readRDS("D:/NMR metabolomics/Raw_Data/ukb_firstreport.rds")

umap_results <- readRDS("D:/NMR metabolomics/Result_Data/pns_umap_lists.rds")

Field_name <- readRDS("D:/NMR metabolomics/Raw_Data/ukb_nmr_colnames.rds")

# Source the function for building manifold
source("ukb_nmr_cox_diseases_build_Manifold.R")


#---------------------------
# Create a list to store all results
#---------------------------
all_manifold_lists <- list()

#---------------------------
# Define c_ids to process
#---------------------------
c_ids <- c(1, 2, 5)

#---------------------------
# Main loop to process each c_id
#---------------------------
for (c_id in c_ids) {
  
  #-----------------------------------------
  # Generate the manifold name (e.g., M1, M2, M5)
  #-----------------------------------------
  m_name <- paste0("M", c_id)
  
  #-----------------------------------------
  # Build the manifold for the current c_id
  #-----------------------------------------
  final_data <- ukb_nmr_cox_diseases_build_Manifold(
    C_id = c_id,
    umap_results = umap_results,
    patient_subset_ID = Participant_ID_Subset_50K,
    ukb_firstreport = ukb_firstreport,
    UKB_All_Data = ukb_all
  )
  
  #-----------------------------------------
  # Read and process the raw data 
  #-----------------------------------------
  # 1) Read in the pre-clustered raw data
  raw_data <- read_csv(sprintf("Fitting_Input_Data/Clustered_raw_C_%d.csv", c_id))
  # 2) Read in the clustering results
  clustering_data <- readRDS(sprintf("Result_Data/M%d_Clustering.rds", c_id))
  # 3) Subset the raw data according to Participant_ID_Subset_50K$data_id, 
  #    remove the first column which might be an index column
  raw_data <- raw_data[Participant_ID_Subset_50K$data_id, -1]
  # 4) Rename the identifier column to match the UKB ID
  raw_data$f.eid <- Participant_ID_Subset_50K$ukb_id
  # 5) Merge raw_data with clustering_data by "f.eid"
  raw_data <- left_join(raw_data, clustering_data, by = "f.eid")
  
  #-----------------------------------------
  # Read and process the fitted data
  #-----------------------------------------
  # 1) Read in the fitted data output from another process
  fitted_data <- read_csv(
    sprintf("Fitting_Output_Data/Final Output/MF_Output_C_%d.csv", c_id), 
    col_names = FALSE
  )
  # 2) Make sure the fitted data column names align with raw_data
  colnames(fitted_data) <- colnames(raw_data)
  # 3) Rename the identifier column for fitted_data as well
  fitted_data$f.eid <- Participant_ID_Subset_50K$ukb_id
  # 4) Merge fitted_data with clustering_data by "f.eid"
  fitted_data <- left_join(fitted_data, clustering_data, by = "f.eid")
  
  #-----------------------------------------
  # Combine and store data in a named list
  #-----------------------------------------
  current_list <- list(
    final = final_data,  # Manifold result 
    raw = raw_data,      # Processed raw data
    fitted = fitted_data # Processed fitted data
  )
  
  #-----------------------------------------
  # Add the current list to the master list
  #-----------------------------------------
  all_manifold_lists[[m_name]] <- current_list
  
  #-----------------------------------------
  # Print progress message
  #-----------------------------------------
  cat(sprintf("Finished processing %s\n", m_name))
}

#---------------------------
# Save all results to an RDS file
#---------------------------
saveRDS(all_manifold_lists, "Result_Data/all_manifold_lists.rds")
cat("All results have been saved to Result_Data/all_manifold_lists.rds\n")

