#' Wrapper function for bootrapping the model fits
#'
#' @param df dataframe
#' @param value_name column name that contains the biomarker value, e.g. "Value"
#' @param time_name column name that contains the time value, e.g. "Time"
#' @param id_name column name that contains the subject id value, e.g. "ID"
#' @param PET_pos_threshold numeric threshold for biomarker positivity
#' @param degree degree of polynomial fit
#' @param num_bootraps number of bootstraps
#' @param bootstrap_percent percentage of IDs to include in each bootstrap iteration
#' @param printIter print the iterations of the bootstrap as a mode of progress
#' @return df_tau dataframe containing ID, predicted value based on model fit, time to positivity estimate
#' @export
#' 
#' 
#' 

bootstrap_get_Time_to_Positivity <- function(df, PET_pos_threshold, id_name, time_name, value_name,
         num_bootstraps = 1000, bootstrap_percent = 0.8, degree = 3, printIter = TRUE) {
  df_res <- list()
  
  Time_Window <- seq(from = -20, to = 20, by = 0.5)
  # List to store bootstrapped results
  df_bs <- vector("list", num_bootstraps)
  
  # Perform bootstrap resampling
  for (i in 1:num_bootstraps) {
    if(printIter == TRUE){
      print(i)
    }
    # Sample IDs with replacement
    sampled_ids <- sample(seq(from = 1, to = length(unique(df[[ id_name]]))), 
                          size = floor(bootstrap_percent * length(unique(df[[id_name]]))), replace = TRUE)
    ID_vec <- unique(df[[id_name]])
    ID_vec <- ID_vec[sampled_ids]
    # Subset the data for the sampled IDs
    df <- data.frame(df)
    sampled_data <- df[df[[id_name]] %in% (ID_vec),] 
    
    # Apply get_Time_to_Positivity on the sampled data
    result <- get_Time_to_Positivity(sampled_data,id_name, time_name, value_name, PET_pos_threshold, degree)
    
    # Store the actual_predicted_val from this iteration
    df_bs[[i]] <- data.frame("val_bs" = result$actual_predicted_val, 
                             "time_bs" = result$Time_to_Positivity)
    interpolated_val <- vector("list", length(Time_Window))
    
    for(j in 1:length(Time_Window)){
      
      interpolated_val[[j]] <- approx( result$Time_to_Positivity, result$actual_predicted_val, as.numeric(Time_Window[j]))$y
    }
    df_res[[i]] <- data.frame("Time_Window" = Time_Window)
    
    df_res[[i]]$interpolated_val <- unlist(interpolated_val)
    
  }
  
  # Convert results to a matrix or dataframe for easier processing
  bootstrap_matrix <- do.call(rbind, df_res)
  mean_result <- data.frame(setDT(bootstrap_matrix)[, median(interpolated_val, na.rm = TRUE), by = Time_Window])
  sd_result <- data.frame(setDT(bootstrap_matrix)[, sd(interpolated_val, na.rm = TRUE), by = Time_Window])
  
  ci_calc <- merge(mean_result, sd_result, by = "Time_Window", suffix = c("_mean", "_sd"))
  
  # Combine into a dataframe of confidence intervals
  ci_df <- data.frame(
    Time_Window = ci_calc$Time_Window,
    Estimate =ci_calc$V1_mean,
    CI_Lower = ci_calc$V1_mean - 1.96 * ci_calc$V1_sd,
    CI_Upper = ci_calc$V1_mean + 1.96 * ci_calc$V1_sd
  )
  # 
  # min_row <- ci_df[ci_df$Estimate < PET_pos_threshold ,] %>%
  #   filter(Estimate == max(Estimate, na.rm = TRUE))
  # 
  # max_row <- ci_df[ci_df$Estimate > PET_pos_threshold ,] %>%
  #   filter(Estimate == min(Estimate, na.rm = TRUE))
  # 
  # adjustment <- approx(ci_df$Estimate, ci_df$Time_Window, PET_pos_threshold)$y
  # ci_df$Time_to_Positivity <- ci_df$Time_Window - as.numeric(adjustment)
  # 
  # return(ci_df[, c("Time_to_Positivity", "Estimate", "CI_Lower", "CI_Upper")])
  
  return(ci_df)
}
