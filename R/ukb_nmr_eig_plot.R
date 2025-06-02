library(ggplot2)

# Function to compute eigenvalues remains the same as provided
compute_eigenvalues <- function(matrix_data) {
  # Convert input to matrix and perform scaling
  matrix_data <- as.matrix(matrix_data)
  scaled_data <- scale(matrix_data, center = TRUE, scale = TRUE)
  # Compute SVD and normalize singular values
  svd_result <- svd(scaled_data)
  singular_values <- svd_result$d
  max_singular_value <- max(singular_values)
  return(singular_values / max_singular_value)
}

# Read raw data and matrix factorization outputs
All_raw <-  read.csv("D:/NMR metabolomics/Fitting_Input_Data/all_biomarker_subset.csv",header = TRUE)[,-1]
C1_RF <- read.csv("D:/NMR metabolomics/Fitting_Output_Data/Final Output/MF_Output_C_1.csv",header = FALSE)
C2_TF <- read.csv("D:/NMR metabolomics/Fitting_Output_Data/Final Output/MF_Output_C_2.csv",header = FALSE)
C3_TF <- read.csv("D:/NMR metabolomics/Fitting_Output_Data/Final Output/MF_Output_C_3.csv",header = FALSE)
C4_RF <- read.csv("D:/NMR metabolomics/Fitting_Output_Data/Final Output/MF_Output_C_4.csv",header = FALSE)
C5_TF <- read.csv("D:/NMR metabolomics/Fitting_Output_Data/Final Output/MF_Output_C_5.csv",header = FALSE)
C6_TF <- read.csv("D:/NMR metabolomics/Fitting_Output_Data/Final Output/MF_Output_C_6.csv",header = FALSE)
C7_TF <- read.csv("D:/NMR metabolomics/Fitting_Output_Data/Final Output/MF_Output_C_7.csv",header = FALSE)

# Calculate eigenvalues for each matrix (take first 10 values)
eigenvalues_all_raw <- compute_eigenvalues(All_raw)[1:10]
eigenvalues_C1_RF <- compute_eigenvalues(C1_RF)[1:10]
eigenvalues_C2_TF <- compute_eigenvalues(C2_TF)[1:10]
eigenvalues_C3_TF <- compute_eigenvalues(C3_TF)[1:10]
eigenvalues_C4_RF <- compute_eigenvalues(C4_RF)[1:10]
eigenvalues_C5_TF <- compute_eigenvalues(C5_TF)[1:10]
eigenvalues_C6_TF <- compute_eigenvalues(C6_TF)[1:10]
eigenvalues_C7_TF <- compute_eigenvalues(C7_TF)[1:10]

# Helper function to create dataframe for plotting
create_eigenvalue_df <- function(singular_values, label) {
  data.frame(EigenValue = singular_values, Matrix = rep(label, length(singular_values)), Index = 1:length(singular_values))
}

# Create individual dataframes for each matrix
df_all_raw <- create_eigenvalue_df(eigenvalues_all_raw, "Raw")
df_C1_RF <- create_eigenvalue_df(eigenvalues_C1_RF, "C1_RF")
df_C2_TF <- create_eigenvalue_df(eigenvalues_C2_TF, "C2_TF")
df_C3_TF <- create_eigenvalue_df(eigenvalues_C3_TF, "C3_TF")
df_C4_RF <- create_eigenvalue_df(eigenvalues_C4_RF, "C4_RF")
df_C5_TF <- create_eigenvalue_df(eigenvalues_C5_TF, "C5_TF")
df_C6_TF <- create_eigenvalue_df(eigenvalues_C6_TF, "C6_TF")
df_C7_TF <- create_eigenvalue_df(eigenvalues_C7_TF, "C7_TF")

#-----------------------------#
#    Subplot Creation Function #
#-----------------------------#
# Function to create individual subplot with consistent styling
create_subplot <- function(data, group_name) {
  ggplot(data, aes(x = Index, y = EigenValue)) +
    geom_line(color = "#4682B4", size = 0.8, alpha = 0.6) +
    geom_point(color = "#00008B", size = 2) +
    labs(
      x = "Index",
      y = "Normalized Eigenvalue",
      title = group_name
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 9, hjust = 0.5),
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 9, color = "black"),
      axis.line = element_line(color = "black", size = 0.5),
      panel.grid = element_blank(),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25)
    ) +
    scale_x_continuous(
      breaks = seq(2, 10, 2)
    )
}

#-----------------------------#
#    Create All Subplots      #
#-----------------------------#
# Prepare data and titles for plotting
plot_data_list <- list(
  df_all_raw, df_C1_RF, df_C2_TF, df_C3_TF, df_C4_RF, 
  df_C5_TF, df_C6_TF, df_C7_TF
)
plot_titles <- c("Raw","M1", "M2", "M3", "M4", 
                 "M5", "M6", "M7")

# Generate all subplots using mapply
plots <- mapply(
  function(data, title) create_subplot(data, title),
  plot_data_list,
  plot_titles,
  SIMPLIFY = FALSE
)

#-----------------------------#
#    Layout Using gridExtra   #
#-----------------------------#
library(gridExtra)
# Arrange plots in a 2x4 grid layout
combined_plot <- grid.arrange(
  grobs = plots,
  ncol = 4,       # 4 columns
  nrow = 2        # 2 rows
)

#-----------------------------#
#    Save Output Figure       #
#-----------------------------#
# Save the combined plot as a high-resolution PDF
ggsave(
  filename = "Plot/eig_plot.pdf",
  plot = combined_plot,
  width = 12,     # Width in inches
  height = 4,     # Height in inches
  units = "in",
  dpi = 300
)
