#' Cross-sectional time to amyloid positivity from raw cohort data
#'
#' This is a one-shot wrapper that builds the cohort t(mu) curve using the
#' existing Budgeon pipeline and then computes time-to-positivity for one or
#' more cross-sectional values.
#'
#' @param df A data.frame with longitudinal cohort data.
#' @param id_name,time_name,value_name Column names (strings) for subject ID, time, and amyloid measure.
#' @param mu_subject Numeric scalar or vector: subject cross-sectional value (same scale as value_name).
#' @param mu_threshold Numeric scalar: positivity cutoff (same scale as value_name).
#' @param degree Integer polynomial degree for slope-vs-level fit (default 3).
#' @return Numeric vector of time-to-positivity (years)
#' @examples
#' \dontrun{
#'   df <- generate_synth_data(150, 25, 1)
#'   t2a_crosssectional(df, "ID", "Time", "Marker_level",
#'                      mu_subject = 0.95, mu_threshold = 1.10)
#' }
#' @export
t2a_crosssectional <- function(df, id_name, time_name, value_name,
                               mu_subject, mu_threshold, degree = 3) {
  # 1) per-ID slope
  ind <- get_individ_lm(df, value_name, time_name, id_name)
  ind <- ind[ind$Slope >= 0, , drop = FALSE]
  if (nrow(ind) < 1L) return(rep(NA_real_, length(mu_subject)))
  
  # 2) eval at midpoint
  eval <- eval_lmm(
    df[df[[id_name]] %in% ind[[id_name]], , drop = FALSE],
    id_name = id_name, time_name = time_name, individual_coeffs = ind
  )
  if (nrow(eval) < 2L) return(rep(NA_real_, length(mu_subject)))
  
  # 3) integrate to get t(mu)
  curve_df <- get_integration_estimates(eval, id_name = id_name, degree = degree)
  
  # 4) cross-sectional T2A+ (vectorized over mu_subject)
  mu_subject <- as.numeric(mu_subject)
  vapply(mu_subject, function(mu)
    time_to_positivity_crosssectional(curve_df, mu, mu_threshold),
    numeric(1)
  )
}
