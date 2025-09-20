#' Generates synthetic data
#'
#' This function creates synthetic datasets to test and demonstrate the Budgeon
#' pipeline. It can generate data under two scenarios controlled by the
#' `baseline_variability` parameter.
#'
#' @param n Number of individuals.
#' @param T Total length of the disease progression timeline.
#' @param l Length of time between observations.
#' @param baseline_variability A character string specifying the baseline scenario.
#'   - `"heterogeneous"` (default): Creates a realistic cohort with varied baseline
#'     biomarker levels. This is the recommended setting for testing the full
#'     pipeline, including reliability filtering.
#'   - `"homogeneous"`: Creates a dataset that replicates the simulation from
#'     Budgeon et al. (2017), where all subjects have similar, low baselines.
#'     Useful for demonstrating the classic model without filtering.
#' @return A synthetic dataset.
#' @export
generate_synth_data <- function (n, T, l, baseline_variability = c("heterogeneous", "homogeneous")){
  baseline_variability <- match.arg(baseline_variability)
  
  b <- runif(n, 10, 16)
  a <- truncnorm::rtruncnorm(n, 0, Inf, 0.4, 0.02) 
  
  if (baseline_variability == "homogeneous") {
    # Homogeneous baselines, replicating the original paper for the classic model
    c <- truncnorm::rtruncnorm(n, 0, Inf, 0.01, 0.005)
  } else {
    # Heterogeneous (varied) baselines for demonstrating the filtering feature
    c <- truncnorm::rtruncnorm(n, 0, Inf, 0.25, 0.15) 
  }
  
  f <- truncnorm::rtruncnorm(n, -Inf, 1, 0.98, 0.005)
  g <- rnorm(n, 20, 1)
  b_star <- b+log(2^(1/g)-1)/a
  tau_all <- seq(0, T, by = l)
  
  all_data <- data.frame()
  for (i in 1: n){
    R_tau <- c[i] + (f[i]-c[i])/(1+exp(-a[i]*(tau_all-b_star[i])))^g[i]
    individual_data<- data.frame(ID = i, Time = tau_all, Marker_level = R_tau)
    all_data <- rbind(all_data, individual_data)
  }
  return(all_data)
}

#' Performs ID-by-ID linear regression
#'
#'
#' @param df dataframe
#' @param value_name column name that contains the biomarker value, e.g. "Value"
#' @param time_name column name that contains the time value, e.g. "Time"
#' @param id_name column name that contains the subject id value, e.g. "ID"

#' @return individual_coeffs dataframe containing slopes and intercepts on an ID-by-ID basis
#' @export

get_individ_lm <- function(df, value_name, time_name, id_name){
  df <- data.frame(df)
  df <- df[!is.na(df[[value_name]]),]
  duplicated_IDs <- unique(df[duplicated(df[[id_name]]), id_name])
 
   
  # df <- function(df){
  #   tryCatch({
      df[df[[id_name]] %in% duplicated_IDs,]
  #   },
  #   error = function(cond){
  #   df[df[[id_name]] %in% duplicated_IDs[[id_name]], ] #alternate depending on package loading situation
  #   })
  # }
  
  subject_col <- sym(id_name)
  
  #Added fix in the event of duplicated Time_from_baselines
  df_filtered <- df %>%
    group_by(!!subject_col) %>%
    filter(n() > 1, n_distinct(!!sym(time_name)) > 1) %>%
    ungroup()
  if (nrow(df_filtered) == 0L) {
    id_empty <- df[0, id_name, drop = FALSE]          
    out <- cbind(id_empty,
                 Intercept = numeric(0),
                 Slope     = numeric(0))
    names(out)[1] <- id_name
    return(out)
  }
  # Create the formula dynamically
  formula <- as.formula(paste(value_name, "~", time_name))
  
  # Run linear regression for each ID separately
  individual_coeffs <- df_filtered %>%
    group_by(!!subject_col) %>%
    group_modify(~ {
      model <- lm(formula, data = .x)
      data.frame(Intercept = coef(model)[1], Slope = coef(model)[2])
    }) %>%
    ungroup()
  
  return(individual_coeffs)
}



#' Calculates mean at midpoint for each individual's trajectory
#'
#'
#' @param df dataframe
#' @param value_name column name that contains the biomarker value, e.g. "Value"
#' @param time_name column name that contains the time value, e.g. "Time"
#' @param id_name column name that contains the subject id value, e.g. "ID"
#' @param individual_coeffs dataframe output from get_individ_lm()

#' @return eval_at_midpoint dataframe containing mean at midpoint of individual trajectory
#' @export
eval_lmm <- function(df, id_name, time_name, individual_coeffs){
  subject_col <- sym(id_name)
  time_col <- sym(time_name)
  df <- df %>%
    group_by(!!subject_col)%>%
    mutate(
      midpoint_time = mean(!!time_col, na.rm = TRUE)
    )
  
  eval_at_midpoint <- df %>%
    left_join(individual_coeffs, by = id_name) %>%
    mutate (
      eval_at_midpoint = Slope * midpoint_time + Intercept
    )
  
  eval_at_midpoint <- eval_at_midpoint %>%
    group_by(!!subject_col)%>%
    select(!!subject_col, midpoint_time, eval_at_midpoint, Slope, Intercept) %>%
    slice(1) 
  return(eval_at_midpoint)
}


#' Fit a non-negative polynomial - helper function
#'
#' @param x time variable
#' @param y biomarker variable
#' @param degree degree of polynomial fit. Budgeon et al suggests 3 is sufficient, but 5 is also good.
#' @return fit polynomial model
#' @export
fit_polynomial <- function(x, y, degree = 3) {
  fit <- lm(y ~ poly(x, degree, raw = TRUE))
  return(fit)
}

#' Fit a non-negative polynomial 
#'
#' @param eval_at_midpoint dataframe containing mean at midpoint of individual trajectory
#' @param id_name column name that contains the subject id value, e.g. "ID"
#' @return to_return dataframe with polynomial fits corresponding with IDs
#' @export

poly_ft <- function(eval_at_midpoint, id_name, degree = 3){
  eval_at_midpoint <- eval_at_midpoint[!is.na(eval_at_midpoint$Slope),]
  x_vals = eval_at_midpoint$eval_at_midpoint
  y_vals = eval_at_midpoint$Slope
  # cubic  
  polynomial_fit <- fit_polynomial(x_vals, y_vals, degree = 3)
  
  # Evaluate the fitted polynomial
  y_fit <- polynomial_fit$fitted.values
  
  to_return <- data.frame(ID = eval_at_midpoint[, id_name],
                          x_vals = x_vals, y_vals=y_vals, y_fit=y_fit)
  colnames(to_return)[1] <- id_name
  return(to_return)
}


#' Shifts model output so that "0" is when the biomarker becomes positive
#' Budgeon et al starts all models at time 0, this shift is done to make the model output consistent with SILA (Betthauser et al, 2022)
#' @param df dataframe containing the IDs, biomarker values and timepoints
#' @param PET_pos_threshold threshold for biomarker positivity
#' @return hor_adj single number used to shift model output time
#' @export
get_horizontal_adjustment <- function(df, PET_pos_threshold){
  if(min(df$actual_predicted_val) < PET_pos_threshold){
    min_row <- df[df$actual_predicted_val < PET_pos_threshold, 
    ] %>% filter(actual_predicted_val == max(actual_predicted_val, na.rm = TRUE))
    max_row <- df[df$actual_predicted_val > PET_pos_threshold, 
    ] %>% filter(actual_predicted_val == min(actual_predicted_val, na.rm = TRUE))
    hor_adj <- ((min_row$tau - max_row$tau) * (PET_pos_threshold - 
                                                 max_row$actual_predicted_val))/(min_row$actual_predicted_val - 
                                                                                   PET_pos_threshold) + max_row$tau} else{
                                                                                     hor_adj <- min(df$actual_predicted_val) - PET_pos_threshold
                                                                                   }
  return(hor_adj)
  
}

#' Integrates per step 3 in Budgeon et al process
#' @param eval_at_midpoint dataframe containing mean at midpoint of individual trajectory
#' @param id_name column name that contains the subject id value, e.g. "ID"
#' @return dataframe with tau and mu_tau (step 4)
#' @export
get_integration_estimates <- function(eval_at_midpoint, id_name, degree = 3){
  polynomial_fit <- poly_ft(eval_at_midpoint = eval_at_midpoint, id_name = id_name, degree = 3)
  polynomial_fit <- polynomial_fit[with(polynomial_fit, order(x_vals)),]
  
  ID <- polynomial_fit[, id_name]
  x_vals <- polynomial_fit$x_vals
  y_vals <- polynomial_fit$y_vals
  y_fit <- polynomial_fit$y_fit
  # Reciprocal
  reciprocal_slope <- 1 / y_fit
  

  
  # integrate the reciprocal slope with respect to fitted mean
  integrated_tau <- cumsum(reciprocal_slope * diff(c(0, x_vals)))
  
  # Find the minimum µ(t) (smallest observed mean)
  min_mu <- min(x_vals)
  mu_tau <- x_vals - min_mu
  
  df_tau <- data.frame(ID = ID, mu_tau = mu_tau, tau = integrated_tau,
                       actual_predicted_val = mu_tau + min_mu)
  
  return(df_tau)
}


#' Wrapper function for entire process
#'
#'
#' @param df dataframe
#' @param value_name column name that contains the biomarker value, e.g. "Value"
#' @param time_name column name that contains the time value, e.g. "Time"
#' @param id_name column name that contains the subject id value, e.g. "ID"
#' @param PET_pos_threshold numeric threshold for biomarker positivity
#' @param degree degree of polynomial fit

#' @return df_tau dataframe containing ID, predicted value based on model fit, time to positivity estimate
#' @export
get_Time_to_Positivity <- function(
    df,
    id_name,
    time_name,
    value_name,
    PET_pos_threshold,
    degree = 3,
    reliability_filter = c("none", "derived_value", "user_value"),
    value_threshold = NULL,
    arc_lower_comp_quantile = 0.99,
    include_negative = TRUE,
    engine = c("cutpointr", "grid"),
    fallback_to_grid = TRUE,
    grid_n = 200L,
    seed = NULL,
    verbose_filter = FALSE
) {
  reliability_filter <- match.arg(reliability_filter)
  engine <- match.arg(engine)
  filter_info <- NULL
  df_use <- df
  
  if (reliability_filter != "none") {
    # 1) Slopes on the full cohort (for ARC labels / training the baseline cut)
    individual_coeffs0 <- get_individ_lm(df, value_name, time_name, id_name)
    
    # 2) Apply chosen filtering mode
    mode_map <- c(derived_value = "derived_value", user_value = "user_value")
    filter_info <- apply_reliability_filter(
      df = df,
      individual_coeffs = individual_coeffs0,
      id_name = id_name,
      time_name = time_name,
      value_name = value_name,
      slope_col = "Slope",
      mode = mode_map[[reliability_filter]],
      value_threshold = value_threshold,
      lower_comp_quantile = arc_lower_comp_quantile,
      include_negative = include_negative,
      use_cutpointr = (engine == "cutpointr"),
      grid_n = grid_n,
      seed = seed,
      verbose = verbose_filter
    )
    
    # 3) Use the filtered data for the downstream Budgeon steps
    df_use <- filter_info$data
  }

  individual_coeffs <- get_individ_lm(df_use, value_name, time_name, id_name)
  individual_coeffs <- individual_coeffs[individual_coeffs$Slope >= 0, ]
  
  eval_at_midpoint <- eval_lmm(
    df_use[df_use[[id_name]] %in% individual_coeffs[[id_name]], ],
    id_name = id_name,
    time_name = time_name,
    individual_coeffs = individual_coeffs
  )
  
  df_tau <- get_integration_estimates(eval_at_midpoint, id_name = id_name, degree = degree)
  
  hor_adj <- get_horizontal_adjustment(df_tau, PET_pos_threshold)
  df_tau$Time_to_Positivity <- df_tau$tau - hor_adj
  
  out <- df_tau[, c("ID", "actual_predicted_val", "Time_to_Positivity")]
  
  if (reliability_filter == "none") {
    return(out)
  } else {
    return(list(result = out, filter_info = filter_info))
  }
}