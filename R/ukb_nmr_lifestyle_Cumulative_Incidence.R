#---------------------------
# Step 1: Load necessary libraries
#---------------------------
library(dbscan)
library(survival)
library(forestplot)
library(ggplot2)
library(ggrepel) 
library(patchwork)

source("ukb_nmr_cox_diseases_build_Manifold.R")

follow_up_end_date <- as.Date("2023-01-01")

#---------------------------
# Step 2: Load data
#---------------------------

ICD10_FieldID <- read.csv("D:/NMR metabolomics/Raw_Data/ICD10_FieldID.csv")

ukb_firstreport <- readRDS("D:/NMR metabolomics/Raw_Data/ukb_firstreport.rds")

umap_results <- readRDS('D:/NMR metabolomics/Result_Data/pns_umap_lists.rds')

Participant_ID_Subset_50K <- read.csv('D:/NMR metabolomics/Raw_Data/Cluster_7C_row_subset.csv')


M1 <- ukb_nmr_cox_diseases_build_Manifold(
  C_id = 1,
  umap_results = umap_results,
  patient_subset_ID = Participant_ID_Subset_50K,
  ukb_firstreport = ukb_firstreport,
  UKB_All_Data = ukb_all
)

M2 <- ukb_nmr_cox_diseases_build_Manifold(
  C_id = 2,
  umap_results = umap_results,
  patient_subset_ID = Participant_ID_Subset_50K,
  ukb_firstreport = ukb_firstreport,
  UKB_All_Data = ukb_all
)

M5 <- ukb_nmr_cox_diseases_build_Manifold(
  C_id = 5,
  umap_results = umap_results,
  patient_subset_ID = Participant_ID_Subset_50K,
  ukb_firstreport = ukb_firstreport,
  UKB_All_Data = ukb_all
)

Mall <- M1 %>%
  # Rename Cluster in M1
  rename(Cluster_M1 = Cluster) %>%
  # Merge with Cluster from M2
  left_join(
    M2 %>% 
      select(f.eid, Cluster_M2 = Cluster),
    by = "f.eid"
  ) %>%
  # Merge with Cluster from M5
  left_join(
    M5 %>% 
      select(f.eid, Cluster_M5 = Cluster),
    by = "f.eid"
  ) %>%
  # Convert clustering results to factors
  mutate(
    Cluster_M1 = factor(Cluster_M1),
    Cluster_M2 = factor(Cluster_M2),
    Cluster_M5 = factor(Cluster_M5)
  )

# Create Cluster_Mix
Mall <- Mall %>%
  mutate(
    # If any cluster is 2, set to 2, otherwise set to 1
    Cluster_Mix = case_when(
      Cluster_M1 == 2 | Cluster_M2 == 2 | Cluster_M5 == 2 ~ 2,
      TRUE ~ 1
    ),
    # Convert to factor
    Cluster_Mix = factor(Cluster_Mix)
  )

# Create a subset where Cluster_Mix = 2
Mall_C2 <- Mall %>%
  filter(Cluster_Mix == 1)

# Load necessary packages
library(dplyr)
library(survival)
library(survminer)

# Find the field_id corresponding to I20
i20_info <- ICD10_FieldID %>%
  filter(ICD10 == "N18")

field_id <- paste0("f.", i20_info$Date.First.Report.Code, ".0.0")

# Data preparation
Mall_C2_analysis <- Mall_C2 %>%
  mutate(
    # Get the disease diagnosis date
    disease_date = as.Date(!!sym(field_id), format = "%Y-%m-%d"),
    # Mark baseline diseases
    is_baseline_disease = !is.na(disease_date) & disease_date <= baseline_date
  ) %>%
  # Remove individuals with baseline diseases
  filter(!is_baseline_disease) %>%
  mutate(
    # Calculate follow-up time: use disease date for those who developed the disease, otherwise use follow-up end date
    event_time = as.numeric(difftime(
      case_when(
        !is.na(disease_date) & disease_date > baseline_date ~ disease_date,
        TRUE ~ follow_up_end_date
      ),
      baseline_date, units = "days"
    )),
    # Event status: 1 for those who developed the disease during follow-up
    event_status = if_else(!is.na(disease_date) & disease_date > baseline_date, 1, 0),
    # Convert lifestyle scores to factors
    sleep_score = factor(sleep_score, levels = c(1, 0), 
                         labels = c("Healthy Sleep", "Unhealthy Sleep")),
    activity_score = factor(activity_score, levels = c(1, 0), 
                            labels = c("Healthy Activity", "Unhealthy Activity")),
    smoking_score = factor(smoking_score, levels = c(1, 0), 
                           labels = c("Non-smoker", "Smoker"))
  )

# Data processing information
cat("Data Processing Information:\n")
cat("Initial sample size:", nrow(Mall_C2), "\n")
cat("Number of individuals with baseline disease:", sum(Mall_C2$is_baseline_disease, na.rm = TRUE), "\n")
cat("Final analysis sample size:", nrow(Mall_C2_analysis), "\n")
cat("Number of cases during follow-up:", sum(Mall_C2_analysis$event_status == 1), "\n")



# Create survival object
surv_obj <- Surv(time = Mall_C2_analysis$event_time, 
                 event = Mall_C2_analysis$event_status)

# Create three survival curves
fit_sleep <- survfit(surv_obj ~ sleep_score, data = Mall_C2_analysis)
fit_activity <- survfit(surv_obj ~ activity_score, data = Mall_C2_analysis)
fit_smoking <- survfit(surv_obj ~ smoking_score, data = Mall_C2_analysis)

# Plot three cumulative incidence curves
p1 <- ggsurvplot(
  fit_sleep,
  data = Mall_C2_analysis,
  title = "Sleep Pattern and Cumulative I20 Incidence",
  fun = "event",
  pval = TRUE,
  conf.int = TRUE,
  palette = c("#2E9FDF", "#E7B800"),
  xlab = "Time (days)",
  ylab = "Cumulative Incidence",
  legend.title = "",
  legend.labs = c("Healthy Sleep", "Unhealthy Sleep"),
  risk.table = FALSE,
  font.title = 10,        # Title font size
  font.x = 8,             # X-axis label font size
  font.y = 8,             # Y-axis label font size
  font.legend = 8,        # Legend font size
  pval.size = 3           # P-value font size
)

p2 <- ggsurvplot(
  fit_activity,
  data = Mall_C2_analysis,
  title = "Physical Activity and Cumulative I20 Incidence",
  fun = "event",
  pval = TRUE,
  conf.int = TRUE,
  palette = c("#2E9FDF", "#E7B800"),
  xlab = "Time (days)",
  ylab = "Cumulative Incidence",
  legend.title = "",
  legend.labs = c("Healthy Activity", "Unhealthy Activity"),
  risk.table = FALSE,
  font.title = 10,
  font.x = 8,
  font.y = 8,
  font.legend = 8,
  pval.size = 3
)

p3 <- ggsurvplot(
  fit_smoking,
  data = Mall_C2_analysis,
  title = "Smoking Status and Cumulative I20 Incidence",
  fun = "event",
  pval = TRUE,
  conf.int = TRUE,
  palette = c("#2E9FDF", "#E7B800"),
  xlab = "Time (days)",
  ylab = "Cumulative Incidence",
  legend.title = "",
  legend.labs = c("Non-smoker", "Smoker"),
  risk.table = FALSE,
  font.title = 10,
  font.x = 8,
  font.y = 8,
  font.legend = 8,
  pval.size = 3
)

# Combine the three plots vertically
p1$plot / p2$plot / p3$plot

# Calculate the final incidence rate for each group
incidence_rates <- Mall_C2_analysis %>%
  summarise(
    # Sleep pattern
    N_Healthy_Sleep = sum(sleep_score == "Healthy Sleep"),
    Cases_Healthy_Sleep = sum(sleep_score == "Healthy Sleep" & event_status == 1),
    Rate_Healthy_Sleep = sprintf("%.2f%%", 
                                 sum(sleep_score == "Healthy Sleep" & event_status == 1) / 
                                   sum(sleep_score == "Healthy Sleep") * 100),
    
    N_Unhealthy_Sleep = sum(sleep_score == "Unhealthy Sleep"),
    Cases_Unhealthy_Sleep = sum(sleep_score == "Unhealthy Sleep" & event_status == 1),
    Rate_Unhealthy_Sleep = sprintf("%.2f%%", 
                                   sum(sleep_score == "Unhealthy Sleep" & event_status == 1) / 
                                     sum(sleep_score == "Unhealthy Sleep") * 100),
    
    # Physical activity
    N_Healthy_Activity = sum(activity_score == "Healthy Activity"),
    Cases_Healthy_Activity = sum(activity_score == "Healthy Activity" & event_status == 1),
    Rate_Healthy_Activity = sprintf("%.2f%%", 
                                    sum(activity_score == "Healthy Activity" & event_status == 1) / 
                                      sum(activity_score == "Healthy Activity") * 100),
    
    N_Unhealthy_Activity = sum(activity_score == "Unhealthy Activity"),
    Cases_Unhealthy_Activity = sum(activity_score == "Unhealthy Activity" & event_status == 1),
    Rate_Unhealthy_Activity = sprintf("%.2f%%", 
                                      sum(activity_score == "Unhealthy Activity" & event_status == 1) / 
                                        sum(activity_score == "Unhealthy Activity") * 100),
    
    # Smoking status
    N_Nonsmoker = sum(smoking_score == "Non-smoker"),
    Cases_Nonsmoker = sum(smoking_score == "Non-smoker" & event_status == 1),
    Rate_Nonsmoker = sprintf("%.2f%%", 
                             sum(smoking_score == "Non-smoker" & event_status == 1) / 
                               sum(smoking_score == "Non-smoker") * 100),
    
    N_Smoker = sum(smoking_score == "Smoker"),
    Cases_Smoker = sum(smoking_score == "Smoker" & event_status == 1),
    Rate_Smoker = sprintf("%.2f%%", 
                          sum(smoking_score == "Smoker" & event_status == 1) / 
                            sum(smoking_score == "Smoker") * 100)
  )

# Print incidence rate information
cat("\nFinal Incidence Rate Statistics:")
cat("\n\nSleep Pattern:")
cat("\nHealthy Sleep Group:", incidence_rates$Rate_Healthy_Sleep, 
    sprintf("(%d/%d)", incidence_rates$Cases_Healthy_Sleep, incidence_rates$N_Healthy_Sleep))
cat("\nUnhealthy Sleep Group:", incidence_rates$Rate_Unhealthy_Sleep,
    sprintf("(%d/%d)", incidence_rates$Cases_Unhealthy_Sleep, incidence_rates$N_Unhealthy_Sleep))

cat("\n\nPhysical Activity:")
cat("\nHealthy Activity Group:", incidence_rates$Rate_Healthy_Activity,
    sprintf("(%d/%d)", incidence_rates$Cases_Healthy_Activity, incidence_rates$N_Healthy_Activity))
cat("\nUnhealthy Activity Group:", incidence_rates$Rate_Unhealthy_Activity,
    sprintf("(%d/%d)", incidence_rates$Cases_Unhealthy_Activity, incidence_rates$N_Unhealthy_Activity))

cat("\n\nSmoking Status:")
cat("\nNon-smoker Group:", incidence_rates$Rate_Nonsmoker,
    sprintf("(%d/%d)", incidence_rates$Cases_Nonsmoker, incidence_rates$N_Nonsmoker))
cat("\nSmoker Group:", incidence_rates$Rate_Smoker,
    sprintf("(%d/%d)", incidence_rates$Cases_Smoker, incidence_rates$N_Smoker))
cat("\n")
