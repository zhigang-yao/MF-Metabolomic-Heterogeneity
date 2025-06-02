###############################################################################
# Set Working Directory
###############################################################################
setwd("D:/NMR metabolomics")

###############################################################################
# Load Required R Packages
###############################################################################
library(uwot)      
library(RANN)      
library(readr)     
library(stats)     
library(ggplot2)   
library(cluster)   
library(cowplot)   
library(dplyr)     
library(dbscan)    
library(tidyr)     
library(gridExtra) 
library(scales)
library(survival)
library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(lubridate)
library(tidyr)
library(patchwork)
library(readxl)
library(readr)
library(writexl)
library(stringr)
library(tidyverse)
library(broom)
source("extract_disease_info.R")

###############################################################################
# Participant Subset and UMAP Results
###############################################################################
Participant_ID_Subset_50K <- read.csv('D:/NMR metabolomics/Raw_Data/Cluster_7C_row_subset.csv')
umap_results <- readRDS('D:/NMR metabolomics/Result_Data/pns_umap_lists.rds')

###############################################################################
# Helper Function: Convert p-value to a Plot Expression
###############################################################################
format_p_expression <- function(pval) {
  if (pval < 2.2e-16) {
    # Very small p-value
    return("italic(P) < 2.2 %*% 10^-16")
  } else {
    # Scientific notation (e.g., "3.1e-05" -> "3.1%*%10^-05")
    val <- formatC(pval, format = "e", digits = 2)
    val <- gsub("e(-?\\d+)", "%*%10^\\1", val)
    return(paste0("italic(P) == ", val))
  }
}

###############################################################################
# 1. Select UMAP Embedding Based on C_id and Perform DBSCAN Clustering
###############################################################################
C_id <- 1
if (C_id %in% c(1, 4, 7)) {
  umap_mainfold_1 <- umap_results$rf_umap_list[[C_id]]
} else {
  umap_mainfold_1 <- umap_results$rec_umap_list[[C_id]]
}

M1 <- data.frame(f.eid = Participant_ID_Subset_50K$ukb_id)
M1 <- cbind(M1, umap_mainfold_1)
colnames(M1)[colnames(M1) == "1"] <- "umap_x"
colnames(M1)[colnames(M1) == "2"] <- "umap_y"

LowDimPoint <- M1[, c("umap_x", "umap_y")]
dbscan_result <- dbscan(LowDimPoint, eps = 1)
M1$Cluster <- dbscan_result$cluster

###############################################################################
# 2. Load and Process Socioeconomic (SES) Data
###############################################################################
ukb_socio <- ukb_all %>%
  select(
    f.eid,
    baseline_date = f.53.0.0, 
    age           = f.21022.0.0,
    sex           = f.31.0.0,
    income        = f.738.0.0,
    education_0   = f.6138.0.0,
    education_1   = f.6138.0.1,
    education_2   = f.6138.0.2,
    education_3   = f.6138.0.3,
    education_4   = f.6138.0.4,
    education_5   = f.6138.0.5,
    career_raw    = f.6142.0.0
  ) %>%
  mutate(
    sex = factor(sex, levels = c(0, 1), labels = c("Female", "Male")),
    baseline_date = as.Date(baseline_date)
  )

ukb_socio$career <- ifelse(ukb_socio$career_raw %in% c(1, 2, 6, 7), 1, 0)

ukb_socio <- ukb_socio %>%
  mutate(
    edu_max = pmax(education_0, education_1, education_2, 
                   education_3, education_4, education_5,
                   na.rm = TRUE),
    edu_max = case_when(
      edu_max == -7 ~ 7,
      edu_max == -3 ~ NA_real_,
      TRUE ~ edu_max
    ),
    SES = case_when(
      # Low
      (career == 0 & edu_max %in% 1:7 & income %in% c(1, 2)) |
        (career == 0 & edu_max == 7 & income %in% c(3, 4, 5)) |
        (career == 1 & income == 1) ~ 1,
      # Mid
      (career == 1 & edu_max %in% 1:6 & income %in% c(2, 3)) |
        (career == 1 & edu_max %in% 3:7 & income == 4) |
        (career == 0 & edu_max %in% 2:6 & income == 3) |
        (career == 0 & edu_max %in% 3:6 & income == 4) ~ 2,
      # High
      TRUE ~ 3
    )
  )

###############################################################################
# 3. Load and Process Physical Measurement Data
###############################################################################
ukb_measure <- ukb_all %>%
  select(
    f.eid,
    baseline_date = f.53.0.0, 
    diastolic_blood_pressure_1 = f.4079.0.0,
    diastolic_blood_pressure_2 = f.4079.0.1,
    systolic_blood_pressure_1  = f.4080.0.0,
    systolic_blood_pressure_2  = f.4080.0.1,
    hand_grip_strength_left    = f.46.0.0,
    hand_grip_strength_right   = f.47.0.0,
    bmi                        = f.21001.0.0,
    c_reactive_protein         = f.30710.0.0
  ) %>%
  mutate(
    diastolic_blood_pressure = (diastolic_blood_pressure_1 + diastolic_blood_pressure_2) / 2,
    systolic_blood_pressure  = (systolic_blood_pressure_1  + systolic_blood_pressure_2)  / 2,
    baseline_date = as.Date(baseline_date)
  )

###############################################################################
# 4. Combine All Data
###############################################################################
all_data <- M1 %>%
  inner_join(ukb_socio,   by = "f.eid") %>%
  inner_join(ukb_measure, by = c("f.eid", "baseline_date")) %>%
  mutate(SES = factor(SES, levels = 1:3, labels = c("Low", "Mid", "High")))

###############################################################################
# 5. Define a Nature-Style Theme
# - P-value font smaller => plot.title = element_text(size=10)
# - x-axis text bigger => axis.text.x = element_text(size=12)
###############################################################################
nature_theme <- theme_minimal() +
  theme(
    plot.background    = element_rect(fill = "white", color = NA),
    panel.grid.major   = element_line(color = "grey90", size = 0.2),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(size = 12, color = "black"),  # Larger bottom axis
    axis.text.y        = element_text(size = 10, color = "black"),  # Keep y-axis text at 10
    axis.title         = element_text(size = 12, color = "black"), 
    plot.title         = element_text(size = 10, face = "bold"),    # Smaller for P-value
    legend.position    = "none",
    plot.margin        = margin(5, 5, 5, 5)
  )

###############################################################################
# 6. Plotting Function #1: Boxplot
###############################################################################
plot_boxplot <- function(data, var_name, var_label, p_value = NULL) {
  q1 <- quantile(data[[var_name]], 0.25, na.rm = TRUE)
  q3 <- quantile(data[[var_name]], 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower <- q1 - 1.5 * iqr
  upper <- q3 + 1.5 * iqr
  
  if (!is.null(p_value)) {
    title_text_expr <- format_p_expression(p_value)
  } else {
    title_text_expr <- "''"
  }
  
  ggplot(data, aes(x = factor(Cluster), y = .data[[var_name]])) +
    geom_boxplot(
      fill = "white",
      color = "#2c3e50", 
      outlier.shape = NA,
      width = 0.5
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 23,
      size = 2,
      fill = "#3498db"
    ) +
    labs(
      title = title_text_expr,
      x = NULL,
      y = var_label
    ) +
    nature_theme +
    scale_x_discrete(labels = as.character(sort(unique(data$Cluster)))) +
    coord_cartesian(
      ylim = c(
        max(lower, min(data[[var_name]], na.rm = TRUE)),
        min(upper, max(data[[var_name]], na.rm = TRUE))
      )
    ) +
    ggtitle(label = parse(text = title_text_expr))  # parse the expression
}

###############################################################################
# 7. Plotting Function #2: Heatmap
###############################################################################
plot_heatmap_transposed <- function(data, var_name, var_label, p_value = NULL) {
  proportions <- data %>%
    group_by(!!sym(var_name), Cluster) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(Cluster) %>%
    mutate(
      proportion = count / sum(count),
      percentage = sprintf("%.1f%%", proportion * 100)
    ) %>%
    ungroup()
  
  if (!is.null(p_value)) {
    title_text_expr <- format_p_expression(p_value)
  } else {
    title_text_expr <- "''"
  }
  
  ggplot(proportions, aes(x = factor(Cluster), y = !!sym(var_name))) +
    geom_tile(aes(fill = proportion), color = "white") +
    geom_text(aes(label = percentage), size = 3, color = "black") +
    scale_fill_gradientn(colors = c("#f7fbff", "#2171b5")) +
    labs(
      title = title_text_expr,
      x = NULL,
      y = var_label
    ) +
    nature_theme +
    scale_x_discrete(labels = as.character(sort(unique(data$Cluster)))) +
    theme(panel.grid = element_blank()) +
    ggtitle(label = parse(text = title_text_expr))
}

###############################################################################
# 8. Generate Plots for Continuous and Categorical Variables
###############################################################################
continuous_vars <- c(
  "age", 
  "systolic_blood_pressure", 
  "diastolic_blood_pressure", 
  "bmi", 
  "hand_grip_strength_left", 
  "hand_grip_strength_right",
  "c_reactive_protein"
)
categorical_vars <- c("sex", "SES")

plots <- list()

# 8.1 Continuous Variables: T-test
for (var in continuous_vars) {
  test <- t.test(as.formula(paste(var, "~ Cluster")), data = all_data)
  pval <- test$p.value
  
  p <- plot_boxplot(
    data      = all_data,
    var_name  = var,
    var_label = var,
    p_value   = pval
  )
  plots[[var]] <- p
}

# 8.2 Categorical Variables: Chi-square
for (var in categorical_vars) {
  test <- chisq.test(table(all_data$Cluster, all_data[[var]]))
  pval <- test$p.value
  
  p <- plot_heatmap_transposed(
    data      = all_data,
    var_name  = var,
    var_label = var,
    p_value   = pval
  )
  plots[[var]] <- p
}

# Update Y-axis labels for the first 7 continuous variables
new_y_labels <- c(
  "Age (years)",
  "SBP (mmHg)",
  "DBP (mmHg)",
  "BMI (kg/m²)",
  "Left Grip (kg)",
  "Right Grip (kg)",
  "CRP (mg/L)"
)

for (i in seq_along(new_y_labels)) {
  plots[[i]] <- plots[[i]] + labs(y = new_y_labels[i])
}

###############################################################################
# 9. Display All Plots
###############################################################################
grid.arrange(grobs = plots, ncol = 9)

###############################################################################
# 10. Create Summary Statistics Table for Continuous Variables
###############################################################################
get_summary_stats <- function(data, var) {
  data %>%
    group_by(Cluster) %>%
    summarise(
      Mean   = mean(!!sym(var), na.rm = TRUE),
      Median = median(!!sym(var), na.rm = TRUE),
      SD     = sd(!!sym(var), na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      Mean   = round(Mean,   2),
      Median = round(Median, 2),
      SD     = round(SD,     2)
    )
}

summary_list <- lapply(continuous_vars, function(var) {
  stats <- get_summary_stats(all_data, var)
  stats$Variable <- var
  stats
})

summary_table <- dplyr::bind_rows(summary_list) %>%
  dplyr::select(Variable, Cluster, Mean, SD, Median) %>%
  dplyr::arrange(Variable, Cluster)

variable_labels <- c(
  "age"                       = "Age (years)",
  "systolic_blood_pressure"   = "SBP (mmHg)",
  "diastolic_blood_pressure"  = "DBP (mmHg)",
  "bmi"                       = "BMI (kg/m²)",
  "hand_grip_strength_left"   = "Left Grip (kg)",
  "hand_grip_strength_right"  = "Right Grip (kg)",
  "c_reactive_protein"        = "CRP (mg/L)"
)

summary_table$Variable <- variable_labels[summary_table$Variable]

library(knitr)
kable(
  summary_table,
  col.names = c("Variable", "Cluster", "Mean", "SD", "Median"),
  caption = "Summary Statistics by Cluster",
  digits = 2,
  format = "html"
)

write.csv(summary_table, "cluster_summary_statistics_M1.csv", row.names = FALSE)
