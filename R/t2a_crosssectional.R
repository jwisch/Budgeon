#' Cross-sectional time to amyloid positivity (with optional reliable-accumulation filtering)
#'
#' Builds the cohort t(mu) curve using the Budgeon pipeline and then computes
#' time-to-positivity for one or more cross-sectional values (mu_subject).
#' Optionally filters the cohort to "reliable accumulators" before building t(mu)
#' using the same three modes as the main wrapper:
#'   - reliability_filter = "none"          : no filtering (original behavior)
#'   - reliability_filter = "derived_value" : ARC -> labels -> max-F1 baseline cut -> filter
#'   - reliability_filter = "user_value"    : filter by a user-supplied baseline cut
#'
#' @param df data.frame with longitudinal cohort data.
#' @param id_name,time_name,value_name column names for subject ID, time, and amyloid measure.
#' @param mu_subject numeric scalar or vector: subject cross-sectional value(s) (same scale as value_name).
#' @param mu_threshold numeric scalar: positivity cutoff (same scale as value_name).
#' @param degree integer polynomial degree for slope-vs-level fit (default 3).
#'
#' @param reliability_filter one of c("none","derived_value","user_value"). Default "none".
#' @param value_threshold numeric or NULL. Required when reliability_filter = "user_value".
#' @param arc_lower_comp_quantile numeric (0,1), quantile for ARC lower component (default 0.99).
#' @param include_negative logical; keep negative slopes for ARC fit (default TRUE).
#' @param engine one of c("cutpointr","grid") for learning the baseline value cut (default "cutpointr").
#' @param fallback_to_grid logical; if TRUE and cutpointr missing, use the internal grid (default TRUE).
#' @param grid_n integer; grid size if grid engine (default 200).
#' @param seed integer or NULL; seed for GMM EM initialization in ARC estimation.
#' @param verbose_filter logical; print messages from the filtering stage (default FALSE).
#'
#' @return Numeric vector of time-to-positivity (years)

#' @export
t2a_crosssectional <- function(
    df, id_name, time_name, value_name,
    mu_subject, mu_threshold, degree = 3,
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
  
  #  (A) filtering on the cohort before building t(mu)
  df_use <- df
  filter_info <- NULL
  
  if (reliability_filter != "none") {
    # Slopes on the full cohort to support ARC labels / baseline-cut learning
    indiv0 <- get_individ_lm(df, value_name, time_name, id_name)
    
    # Map to apply_reliability_filter() modes
    mode_map <- c(derived_value = "derived_value", user_value = "user_value")
    
    filter_info <- apply_reliability_filter(
      df = df,
      individual_coeffs = indiv0,
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
    df_use <- filter_info$data
  }
  
  # (B) Build the cohort t(mu) curve from the (possibly filtered) data
  ind <- get_individ_lm(df_use, value_name, time_name, id_name)
  ind <- ind[ind$Slope >= 0, , drop = FALSE]
  if (nrow(ind) < 1L) {
    res <- rep(NA_real_, length(mu_subject))
    return(if (reliability_filter == "none") res else list(result = res, curve_df = NULL, filter_info = filter_info))
  }
  
  eval <- eval_lmm(
    df_use[df_use[[id_name]] %in% ind[[id_name]], , drop = FALSE],
    id_name = id_name, time_name = time_name, individual_coeffs = ind
  )
  if (nrow(eval) < 2L) {
    res <- rep(NA_real_, length(mu_subject))
    return(if (reliability_filter == "none") res else list(result = res, curve_df = NULL, filter_info = filter_info))
  }
  
  curve_df <- get_integration_estimates(eval, id_name = id_name, degree = degree)
  
  # (C) Cross-sectional T2A+ for each mu_subject
  mu_subject <- as.numeric(mu_subject)
  res <- vapply(
    mu_subject,
    function(mu) time_to_positivity_crosssectional(curve_df, mu, mu_threshold),
    numeric(1)
  )
  
  if (reliability_filter == "none") {
    res
  } else {
    list(result = res, curve_df = curve_df, filter_info = filter_info)
  }
}
