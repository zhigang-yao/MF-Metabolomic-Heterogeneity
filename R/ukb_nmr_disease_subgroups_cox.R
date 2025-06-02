#---------------------------
# Step 1: Load required libraries
#---------------------------
library(dbscan)
library(survival)
library(forestplot)
library(ggplot2)
library(ggrepel) 
source("Final_pnas_nmr_cox_diseases_build_M.R")

#---------------------------
# Step 2: Load data
#---------------------------
UKB_All_Data <- readRDS("D:/NMR metabolomics/UKB_All_Data.rds")

patient_subset_ID <- read.csv('Cluster_7C_row_subset.csv')

ICD10_FieldID <- read.csv("ICD10_FieldID.csv")

ukb_firstreport <- readRDS("ukb_firstreport.rds")

umap_results <- readRDS('pns_umap_lists.rds')


M1 <- Final_pnas_nmr_cox_diseases_build_M(
  C_id = 5,
  umap_results = umap_results,
  patient_subset_ID = patient_subset_ID,
  ukb_firstreport = ukb_firstreport,
  UKB_All_Data = UKB_All_Data
)

follow_up_end_date <- as.Date("2023-1-1")


# Modified code
results <- data.frame(
  Disease = character(),
  Disease_Name = character(),
  HR = numeric(),
  Lower_CI = numeric(),
  Upper_CI = numeric(),
  P_Value = numeric(),
  Recall = numeric(),
  stringsAsFactors = FALSE
)

high_risk_cluster <- 2

# Create progress bar
total_diseases <- nrow(ICD10_FieldID)
pb <- txtProgressBar(min = 0, max = total_diseases, style = 3)

# Iterate over each row of ICD10_FieldID
for(i in 1:nrow(ICD10_FieldID)) {
  setTxtProgressBar(pb, i)
  
  disease_code <- ICD10_FieldID[i, 1]  
  disease_name <- ICD10_FieldID[i, 3]
  field_id <- paste0("f.", ICD10_FieldID[i, 2], ".0.0")
  
  if (field_id %in% colnames(M1)) {
    M1[[field_id]] <- as.Date(M1[[field_id]], format = "%Y-%m-%d")
    
    event_time <- pmin(
      difftime(as.Date(ifelse(is.na(M1[[field_id]]), follow_up_end_date, M1[[field_id]])), 
               M1$baseline_date, units = "days")
    )
    
    event_status <- ifelse(is.na(M1[[field_id]]), 0, 1)
    
    # Construct survival object
    surv_obj <- Surv(time = as.numeric(event_time), event = event_status)
    cox_model <- coxph(surv_obj ~ Cluster, data = M1)
    summary_model <- summary(cox_model)
    
    HR <- summary_model$coef[1, 2]
    Lower_CI <- summary_model$conf.int[1, 3]
    Upper_CI <- summary_model$conf.int[1, 4]
    P_Value <- summary_model$coef[1, 5]
    
    total_disease_cases <- sum(event_status == 1)
    cases_in_cluster <- sum(event_status == 1 & M1$Cluster == high_risk_cluster)
    recall <- ifelse(total_disease_cases > 0, cases_in_cluster / total_disease_cases, 0)
    
    # Create result for current disease
    current_result <- data.frame(
      Disease = disease_code,
      Disease_Name = disease_name,
      HR = round(HR, 3),
      Lower_CI = round(Lower_CI, 3),
      Upper_CI = round(Upper_CI, 3),
      P_Value = round(P_Value, 4),
      Recall = round(recall, 3),
      Total_Cases = total_disease_cases
    )
    
    # Append result to total results
    results <- rbind(results, current_result)
    
    # Output current disease result immediately
    cat("\nAnalysis complete:", disease_code, "-", disease_name, "\n")
    cat("HR:", round(HR, 3), 
        sprintf("(95%% CI: %.3f-%.3f)", Lower_CI, Upper_CI),
        "| P-value:", format(P_Value, scientific = TRUE, digits = 4),
        "| Recall:", round(recall, 3),
        "| Total cases:", total_disease_cases, "\n")
    cat("----------------------------------------\n")
    
    # Update CSV file after each disease
    write.csv(results, "disease_analysis_results.csv", row.names = FALSE)
  }
}

close(pb)


# Final overall statistics output
cat("\n=== Analysis complete! ===\n")
cat("\nBasic statistics:\n")
cat("Total diseases analyzed:", nrow(results), "\n")
cat("Diseases with HR > 1:", sum(results$HR > 1, na.rm = TRUE), "\n")
cat("Diseases with HR < 1:", sum(results$HR < 1, na.rm = TRUE), "\n")
cat("Diseases with P < 0.05:", sum(results$P_Value < 0.05, na.rm = TRUE), "\n")
cat("Diseases with Recall > 0.5:", sum(results$Recall > 0.5, na.rm = TRUE), "\n")

cat("\nCombined statistics:\n")
# Combination of HR and P-value
cat("Diseases with HR > 1 and P < 0.05:", 
    sum(results$HR > 1 & results$P_Value < 0.05, na.rm = TRUE), "\n")
cat("Diseases with HR < 1 and P < 0.05:", 
    sum(results$HR < 1 & results$P_Value < 0.05, na.rm = TRUE), "\n")

# Combination of HR and Recall
cat("Diseases with HR > 1 and Recall > 0.5:", 
    sum(results$HR > 1 & results$Recall > 0.5, na.rm = TRUE), "\n")
cat("Diseases with HR < 1 and Recall > 0.5:",
    sum(results$HR < 1 & results$Recall > 0.5, na.rm = TRUE), "\n")

# Combination of P-value and Recall
cat("Diseases with P < 0.05 and Recall > 0.5:", 
    sum(results$P_Value < 0.05 & results$Recall > 0.5, na.rm = TRUE), "\n")

cat("\nThree-indicator combination statistics:\n")
# Combination of all three indicators
cat("Diseases with HR > 1, P < 0.05, and Recall > 0.5:", 
    sum(results$HR > 1 & results$P_Value < 0.05 & results$Recall > 0.5, na.rm = TRUE), "\n")
cat("Diseases with HR < 1, P < 0.05, and Recall > 0.5:", 
    sum(results$HR < 1 & results$P_Value < 0.05 & results$Recall > 0.5, na.rm = TRUE), "\n")


cat("\nPercentage statistics:\n")
total <- nrow(results)
# Calculate the percentage of key combinations
cat("Proportion of significant results (P < 0.05):", 
    round(sum(results$P_Value < 0.05, na.rm = TRUE) / total * 100, 2), "%\n")
cat("Proportion of high-risk significant results (HR > 1 and P < 0.05):", 
    round(sum(results$HR > 1 & results$P_Value < 0.05, na.rm = TRUE) / total * 100, 2), "%\n")
cat("Proportion of protective significant results (HR < 1 and P < 0.05):", 
    round(sum(results$HR < 1 & results$P_Value < 0.05, na.rm = TRUE) / total * 100, 2), "%\n")


# Subsequent significant filtering, plotting, etc., remain unchanged
significant_results <- results[results$P_Value < 0.05, ]
top_30_diseases <- significant_results[order(significant_results$HR, decreasing = TRUE), ][1:30, ]

library(ggplot2)

# Assuming top_30_diseases already contains the necessary columns, sorted by HR
top_30_diseases <- top_30_diseases[order(top_30_diseases$HR, decreasing = TRUE), ]

# Log10 transformation of HR and CI values
top_30_diseases$HR_log <- log10(top_30_diseases$HR)
top_30_diseases$Lower_CI_log <- log10(top_30_diseases$Lower_CI)
top_30_diseases$Upper_CI_log <- log10(top_30_diseases$Upper_CI)

# Calculate the range of log-transformed values
min_HR_log <- min(top_30_diseases$Lower_CI_log, na.rm = TRUE)
max_HR_log <- max(top_30_diseases$Upper_CI_log, na.rm = TRUE)

# Scale log-transformed values to 0-1 range
top_30_diseases$HR_scaled <- (top_30_diseases$HR_log - min_HR_log) / (max_HR_log - min_HR_log)
top_30_diseases$Lower_CI_scaled <- (top_30_diseases$Lower_CI_log - min_HR_log) / (max_HR_log - min_HR_log)
top_30_diseases$Upper_CI_scaled <- (top_30_diseases$Upper_CI_log - min_HR_log) / (max_HR_log - min_HR_log)

# Add significance markers
top_30_diseases$sig_symbols <- ifelse(top_30_diseases$P_Value < 0.001, "***",
                                      ifelse(top_30_diseases$P_Value < 0.01, "**",
                                             ifelse(top_30_diseases$P_Value < 0.05, "*", "")))

# Plotting
ggplot(top_30_diseases, aes(x = reorder(Disease, HR))) +
  geom_errorbar(aes(ymin = Lower_CI_scaled, ymax = Upper_CI_scaled, color = "HR"), width = 0.2) +
  geom_point(aes(y = HR_scaled, color = "HR", shape = "HR"), size = 3) +
  geom_point(aes(y = Recall, color = "Recall", shape = "Recall"), size = 3) +
  # Add significance markers
  geom_text(aes(y = Upper_CI_scaled + 0.1, label = sig_symbols), 
            color = "black", size = 4, show.legend = FALSE) +
  scale_y_continuous(
    name = "Hazard Ratio (HR)",
    limits = c(0, 1.1),  # Slightly expand y-axis to accommodate significance markers
    labels = function(x) format(10^(x * (max_HR_log - min_HR_log) + min_HR_log), digits = 2, scientific = FALSE),
    sec.axis = sec_axis(~., name = "Recall")
  ) +
  xlab("Disease") +
  scale_color_manual(
    values = c("HR" = "black", "Recall" = "forestgreen"),
    name = NULL
  ) +
  scale_shape_manual(
    values = c("HR" = 19, "Recall" = 17),
    name = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.title.y.right = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  )
