ukb_nmr_cox_diseases_build_Manifold <- function(C_id, umap_results, patient_subset_ID, ukb_firstreport, UKB_All_Data) {
  # Step 1: Get UMAP results based on C_id
  umap_mainfold_1 <- if (C_id %in% c(1,4)) {
    umap_results$rf_umap_list[[C_id]]
  } else {
    umap_results$rec_umap_list[[C_id]]
  }
  
  # Step 2: Create initial M1 dataframe and merge with UMAP results
  M1 <- data.frame(f.eid = patient_subset_ID$ukb_id)
  ukb_firstreport_subset <- ukb_firstreport[ukb_firstreport$f.eid %in% M1$f.eid, ]
  M1 <- cbind(M1, umap_mainfold_1)
  colnames(M1)[colnames(M1) == "1"] <- "umap_x"
  colnames(M1)[colnames(M1) == "2"] <- "umap_y"
  M1 <- merge(M1, ukb_firstreport_subset, by = "f.eid")
  
  # Step 3: Perform DBSCAN clustering
  LowDimPoint <- M1[, c("umap_x", "umap_y")]
  dbscan_result <- dbscan(LowDimPoint, eps = 1)
  M1$Cluster <- dbscan_result$cluster
  
  # Step 4: Select baseline measurement data and create new variables
  Baseline_Measurement <- UKB_All_Data %>%
    select(
      f.eid,
      baseline_date = f.53.0.0,
      age = f.21022.0.0,
      sex = f.31.0.0,
      BMI = f.21001.0.0,
      income = f.738.0.0,
      education_0 = f.6138.0.0,
      education_1 = f.6138.0.1,
      education_2 = f.6138.0.2,
      education_3 = f.6138.0.3,
      education_4 = f.6138.0.4,
      education_5 = f.6138.0.5,
      career_raw = f.6142.0.0,
      smoking_1 = f.1239.0.0,
      smoking_2 = f.1249.0.0,
      smoking_3 = f.2644.0.0,
      activity_1 = f.22038.0.0,
      activity_2 = f.22039.0.0,
      sleep_raw = f.1160.0.0,
      salt_added = f.1478.0.0
    ) %>%
    mutate(
      sex = factor(sex, levels = c(0, 1), labels = c("Female", "Male")),
      
      # Create basic categorical variables
      smoking = case_when(
        (smoking_1 == 0 & smoking_2 == 4) |
          (smoking_1 == 0 & smoking_2 %in% c(2,3) & smoking_3 == 0) ~ "non-smoker",
        TRUE ~ "smoker"
      ),
      
      activity = case_when(
        (activity_1 + activity_2) <= 600 ~ "low-activity",
        (activity_1 + activity_2) > 600 & (activity_1 + activity_2) <= 1200 ~ "mid-activity",
        (activity_1 + activity_2) > 1200 ~ "high-activity"
      ),
      
      sleep = case_when(
        sleep_raw < 7 ~ "short-sleep",
        sleep_raw >= 7 & sleep_raw <= 9 ~ "normal-sleep",
        sleep_raw > 9 ~ "long-sleep"
      ),
      
      # Calculate health scores
      sleep_score = case_when(
        sleep == "normal-sleep" ~ 1,
        TRUE ~ 0
      ),
      
      activity_score = case_when(
        activity %in% c("mid-activity", "high-activity") ~ 1,
        TRUE ~ 0
      ),
      
      smoking_score = case_when(
        smoking == "non-smoker" ~ 1,
        TRUE ~ 0
      ),
      
      salt_score = case_when(
        salt_added %in% c(1, 2) ~ 1,
        TRUE ~ 0
      ),
      
      # Calculate the total health score
      health_total_score = sleep_score + activity_score + smoking_score,
      
      # Classify health status
      health_category = case_when(
        health_total_score <= 0 ~ "Unhealthy",
        health_total_score <= 2 ~ "Moderate",
        health_total_score == 3 ~ "Healthy"
      ),
      
      # Convert to ordinal factor
      health_category = factor(
        health_category, 
        levels = c("Unhealthy", "Moderate", "Healthy"),
        ordered = TRUE
      ),
      
      # Continue with original education and SES calculations
      edu_max = pmax(education_0, education_1, education_2, education_3, education_4, education_5, na.rm = TRUE),
      edu_max = case_when(
        edu_max == -7 ~ 7,
        edu_max == -3 ~ NA_real_,
        TRUE ~ edu_max
      ),
      SES = case_when(
        (career_raw %in% c(1, 2, 6, 7) & edu_max %in% 1:7 & income %in% c(1, 2)) |
          (career_raw %in% c(1, 2, 6, 7) & edu_max == 7 & income %in% c(3, 4, 5)) ~ 1,
        (career_raw %in% c(1, 2, 6, 7) & edu_max %in% 1:6 & income %in% c(2, 3)) |
          (career_raw %in% c(1, 2, 6, 7) & edu_max %in% 3:7 & income == 4) ~ 2,
        TRUE ~ 3
      )
    )
  
  # Step 5: Merge all data
  M1 <- M1 %>%
    inner_join(Baseline_Measurement, by = "f.eid") %>%
    mutate(
      SES = factor(SES, levels = 1:3, labels = c("Low", "Mid", "High")),
      activity = factor(activity, levels = c("low-activity", "mid-activity", "high-activity")),
      sleep = factor(sleep, levels = c("short-sleep", "normal-sleep", "long-sleep")),
      smoking = factor(smoking, levels = c("non-smoker", "smoker"))
    )
  
  return(M1)
}
