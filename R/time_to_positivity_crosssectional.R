#' Cross-sectional time-to-positivity
#' @param fit_df Data frame from get_integration_estimates or get_Time_to_Positivity
#' @param mu_subject Numeric: subject's amyloid value
#' @param mu_threshold Numeric: positivity cutoff
#' @return Numeric: estimated years until positivity
#' @export
time_to_positivity_crosssectional <- function(fit_df, mu_subject, mu_threshold) {
  fit_df <- fit_df[order(fit_df$actual_predicted_val), ]
  
  # interpolation function t(mu)
  f_t <- stats::approxfun(fit_df$actual_predicted_val, fit_df$tau)
  
  t_subject <- f_t(mu_subject)
  t_threshold <- f_t(mu_threshold)
  
  t_threshold - t_subject
}
