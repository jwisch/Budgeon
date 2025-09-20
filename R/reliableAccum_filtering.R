#' Estimates a reliable accumulation threshold using 2-component GMM on 
#' the distribution of ARC.
#' 
#' @param individual_coeffs data.frame/tibble returned by get_individ_lm(); must contain {id_name, slope_col}
#' @param id_name character, name of ID column (default "ID")
#' @param slope_col character, name of slope column (default "Slope")
#' @param lower_comp_quantile numeric in (0,1), quantile of the lower component (default 0.99)
#' @param include_negative logical, whether to keep negative slopes in the fit (default TRUE; recommended)
#' @param seed integer or NULL, reproducibility seed for EM initialization (default NULL)
#' @param verbose logical, emit informative messages (default FALSE)
#' @return list with elements:
#' \describe{
#'   \item{arc_threshold}{numeric; estimated threshold in ARC units (e.g., SUVR/year, CL/year)}
#'   \item{labels}{data.frame with {ID, Slope, is_reliable}}
#'   \item{mixture}{list with mixture parameters (means, sds, weights, method, quantile)}
#'   \item{n_used}{number of subjects used}
#'   \item{n_reliable}{number classified as reliable}
#'   \item{pct_reliable}{percentage classified as reliable}
#'   \item{warnings}{character vector of warnings (possibly empty)}
#' }
#' @export
compute_arc_threshold <- function(
    individual_coeffs,
    id_name = "ID",
    slope_col = "Slope",
    lower_comp_quantile = 0.99,
    include_negative = TRUE,
    seed = 123,
    verbose = FALSE
) {
  warn_vec <- character(0)
  
  stopifnot(is.data.frame(individual_coeffs))
  if (!all(c(id_name, slope_col) %in% names(individual_coeffs))) {
    stop("the input dataframe must contain columns: ", id_name, " and ", slope_col)
  }
  
  # Prepare slopes
  slopes <- individual_coeffs[[slope_col]]
  slopes <- slopes[is.finite(slopes)]
  if (!include_negative) {
    slopes <- pmax(slopes, 0)
    if (verbose) message("Negative slopes were clipped to 0 (include_negative=FALSE).")
  }
  slopes <- stats::na.omit(slopes)
  n_used <- length(slopes)
  
  if (n_used < 5L) {
    stop("Not enough subjects to estimate an ARC threshold (n < 5).")
  }
  
  # Empirical fallback (nonparametric): 99th pct within "low-ARC band"
  empirical_threshold <- function(x, q = lower_comp_quantile) {
    band_cut <- stats::quantile(x, 0.60, na.rm = TRUE, names = FALSE)
    low_band <- x[x <= band_cut]
    if (length(low_band) < 5L) {
      warn_vec <<- c(warn_vec, "Empirical fallback: 'low-ARC band' too small; using overall quantile.")
      return(stats::quantile(x, q, na.rm = TRUE, names = FALSE))
    }
    stats::quantile(low_band, q, na.rm = TRUE, names = FALSE)
  }
  
  # If very small sample, skip GMM
  if (n_used < 30L) {
    warn_vec <- c(warn_vec, sprintf("n=%d < 30; using empirical fallback (no GMM).", n_used))
    thres <- empirical_threshold(slopes, q = lower_comp_quantile)
    labels <- individual_coeffs[, c(id_name, slope_col), drop = FALSE]
    labels$is_reliable <- is.finite(labels[[slope_col]]) & (labels[[slope_col]] >= thres)
    return(list(
      arc_threshold = as.numeric(thres),
      labels = labels,
      mixture = list(
        means = NA_real_, sds = NA_real_, weights = NA_real_,
        method = "empirical", lower_comp_quantile = lower_comp_quantile
      ),
      n_used = n_used,
      n_reliable = sum(labels$is_reliable, na.rm = TRUE),
      pct_reliable = 100 * mean(labels$is_reliable, na.rm = TRUE),
      warnings = unique(warn_vec)
    ))
  }
  
  # Try GMM; otherwise fallback
  used_method <- "empirical"
  mix_means <- mix_sds <- mix_weights <- NA_real_
  thres <- NA_real_
  
  if (requireNamespace("mixtools", quietly = TRUE)) {
    if (!is.null(seed)) set.seed(seed)
    fit <- try(mixtools::normalmixEM(slopes, k = 2, maxit = 10000, epsilon = 1e-08), silent = TRUE)
    
    if (!inherits(fit, "try-error") && isTRUE(fit$conv)) {
      idx_low <- which.min(fit$mu)
      mu_low <- fit$mu[idx_low]
      sd_low <- fit$sigma[idx_low]
      lam <- fit$lambda
      
      # Separation diagnostic (warn if weak)
      pooled_sd <- sqrt(mean(fit$sigma^2))
      sep <- abs(diff(sort(fit$mu))) / pooled_sd
      if (sep < 0.25) {
        warn_vec <- c(warn_vec, "GMM: weak separation (|Δμ|/pooled_sd < 0.25). ")
      }
      
      # 99th percentile of lower component
      thres <- stats::qnorm(p = lower_comp_quantile, mean = mu_low, sd = sd_low)
      mix_means <- fit$mu
      mix_sds <- fit$sigma
      mix_weights <- lam
      used_method <- "gmm"
    } else {
      warn_vec <- c(warn_vec, "GMM fit failed or did not converge; using empirical fallback.")
      thres <- empirical_threshold(slopes, q = lower_comp_quantile)
    }
  } else {
    warn_vec <- c(warn_vec, "Package 'mixtools' not available; using empirical fallback.")
    thres <- empirical_threshold(slopes, q = lower_comp_quantile)
  }
  
  labels <- individual_coeffs[, c(id_name, slope_col), drop = FALSE]
  labels$is_reliable <- is.finite(labels[[slope_col]]) & (labels[[slope_col]] >= thres)
  
  res <- list(
    arc_threshold = as.numeric(thres),
    labels = labels,
    mixture = list(
      means = as.numeric(mix_means),
      sds = as.numeric(mix_sds),
      weights = as.numeric(mix_weights),
      method = used_method,
      lower_comp_quantile = lower_comp_quantile
    ),
    n_used = n_used,
    n_reliable = sum(labels$is_reliable, na.rm = TRUE),
    pct_reliable = 100 * mean(labels$is_reliable, na.rm = TRUE),
    warnings = unique(warn_vec)
  )
  if (verbose) {
    message(sprintf(
      "ARC threshold = %.6f (method: %s). Reliable: %d / %d (%.1f%%).",
      res$arc_threshold, res$mixture$method, res$n_reliable, n_used, res$pct_reliable
    ))
    if (length(res$warnings)) message(paste("Note:", paste(res$warnings, collapse = " | ")))
  }
  res
}

#' Compute a baseline-value cutpoint (max F1) using earliest-visit baseline
#'
#' @param df data.frame with {id_name, time_name, value_name}
#' @param arc_labels data.frame with {id_name, is_reliable}
#' @param id_name,time_name,value_name column names in df (defaults: "ID","Time","Marker_level")
#' @param engine one of c("cutpointr","grid"). Default "cutpointr".
#' @param fallback_to_grid logical; if TRUE and engine="cutpointr" but the package is missing,
#'        fall back to internal grid search. Default FALSE (stop instead).
#' @param install_if_missing logical; if TRUE and engine="cutpointr" but package is missing,
#'        attempt utils::install.packages("cutpointr") then continue. Default TRUE
#' @param grid_n integer, grid size if using the internal grid engine (default 200)
#' @param verbose logical
#'
#' @return list: value_threshold, summary (F1, precision, recall, specificity, TP/FP/FN/TN),
#'         baseline_table, warnings
#' @export
compute_value_threshold_via_f1 <- function(
    df,
    arc_labels,
    id_name = "ID",
    time_name = "Time",
    value_name = "Marker_level",
    engine = c("cutpointr", "grid"),
    fallback_to_grid = TRUE,
    install_if_missing = FALSE,
    grid_n = 200L,
    verbose = FALSE
) {
  engine <- match.arg(engine)
  warn_vec <- character(0)
  
  stopifnot(is.data.frame(df), is.data.frame(arc_labels))
  if (!all(c(id_name, time_name, value_name) %in% names(df))) {
    stop("`df` must contain columns: ", paste(c(id_name, time_name, value_name), collapse = ", "))
  }
  if (!all(c(id_name, "is_reliable") %in% names(arc_labels))) {
    stop("`arc_labels` must contain columns: ", id_name, ", is_reliable")
  }
  
  #baseline marker level
  df_local <- df[, c(id_name, time_name, value_name)]
  oo <- order(df_local[[id_name]], df_local[[time_name]])
  df_local <- df_local[oo, , drop = FALSE]
  first_rows <- !duplicated(df_local[[id_name]])  # first row per ID after sorting = earliest time
  baseline <- df_local[first_rows, c(id_name, value_name), drop = FALSE]
  names(baseline)[2] <- "baseline_value"
  
  fit_tab <- merge(baseline, arc_labels[, c(id_name, "is_reliable")], by = id_name, all.x = TRUE)
  fit_tab <- fit_tab[is.finite(fit_tab$baseline_value) & !is.na(fit_tab$is_reliable), , drop = FALSE]
  if (nrow(fit_tab) < 10L) {
    warn_vec <- c(warn_vec, sprintf("Only %d labeled IDs available for baseline cutpoint; results may be unstable.", nrow(fit_tab)))
  }
  
  # Metrics 
  metrics_at <- function(thr) {
    pred <- fit_tab$baseline_value >= thr
    TP <- sum(pred & fit_tab$is_reliable, na.rm = TRUE)
    FP <- sum(pred & !fit_tab$is_reliable, na.rm = TRUE)
    FN <- sum(!pred & fit_tab$is_reliable, na.rm = TRUE)
    TN <- sum(!pred & !fit_tab$is_reliable, na.rm = TRUE)
    precision <- if ((TP + FP) > 0) TP / (TP + FP) else 0
    recall    <- if ((TP + FN) > 0) TP / (TP + FN) else 0
    specificity <- if ((TN + FP) > 0) TN / (TN + FP) else 0
    F1 <- if ((precision + recall) > 0) (2 * precision * recall) / (precision + recall) else 0
    list(F1 = F1, precision = precision, recall = recall,
         sensitivity = recall, specificity = specificity,
         TP = TP, FP = FP, FN = FN, TN = TN)
  }
  
  if (engine == "cutpointr") {
    has_cp <- requireNamespace("cutpointr", quietly = TRUE)
    if (!has_cp && install_if_missing) {
      if (verbose) message("Installing 'cutpointr' ...")
      utils::install.packages("cutpointr")
      has_cp <- requireNamespace("cutpointr", quietly = TRUE)
    }
    if (!has_cp && !fallback_to_grid) {
      stop("Package 'cutpointr' is required for engine='cutpointr'. ",
           "Install it via install.packages('cutpointr') ",
           "or call with install_if_missing=TRUE, or set engine='grid'.")
    }
    if (has_cp) {
      if (verbose) message("Using 'cutpointr' to maximize F1 for baseline cutpoint (earliest visit).")
      d <- fit_tab
      cp <- cutpointr::cutpointr(
        data   = d,
        x      = baseline_value,
        class  = is_reliable,
        method = cutpointr::maximize_metric,
        metric = cutpointr::F1_score
      )
      thr <- as.numeric(cp$optimal_cutpoint)
      met <- metrics_at(thr)
    } else {
      warn_vec <- c(warn_vec, "Package 'cutpointr' not available; falling back to internal grid search.")
      engine <- "grid"
    }
  }
  
  if (engine == "grid") {
    vals <- sort(unique(fit_tab$baseline_value))
    cand <- if (length(vals) <= grid_n) vals else as.numeric(stats::quantile(vals, probs = seq(0, 1, length.out = grid_n)))
    best <- list(F1 = -Inf, thr = NA_real_, met = NULL)
    for (t in cand) {
      m <- metrics_at(t)
      if (m$F1 > best$F1 || (isTRUE(all.equal(m$F1, best$F1)) && t < best$thr)) {
        best <- list(F1 = m$F1, thr = t, met = m)
      }
    }
    thr <- best$thr
    met <- best$met
    if (verbose) message("Used internal grid search to maximize F1 for baseline cutpoint.")
  }
  
  if (verbose) {
    message(sprintf(
      "Baseline value cutpoint (max F1 at earliest visit): %.6f | F1=%.3f, P=%.3f, R=%.3f, Spec=%.3f",
      thr, met$F1, met$precision, met$recall, met$specificity
    ))
  }
  
  list(
    value_threshold = thr,
    summary = met,
    baseline_table = fit_tab,
    warnings = unique(warn_vec)
  )
}


#' Apply reliability filtering with exactly three options:
#'  1) "derived_value": learn baseline-value cutpoint via ARC labels + max-F1, then filter
#'  2) "user_value":    use a user-supplied baseline-value cutpoint, then filter
#'  3) "none":          return the data unchanged
#'
#' @param df long data.frame with columns {id_name, time_name, value_name}
#' @param individual_coeffs data.frame from get_individ_lm(); must include {id_name, slope_col}
#' @param id_name,time_name,value_name column names in df (defaults "ID","Time","Marker_level")
#' @param slope_col column name for the slope in individual_coeffs (default "Slope")
#' @param mode one of c("derived_value","user_value","none")
#' @param value_threshold numeric or NULL.
#'   - required when mode = "user_value" (baseline cut to apply, same units as value_name)
#'   - ignored when mode = "derived_value" or "none"
#' @param lower_comp_quantile numeric in (0,1) for ARC threshold estimation (default 0.99)
#' @param include_negative logical for ARC estimation (default TRUE)
#' @param use_cutpointr logical, prefer cutpointr for F1 (default TRUE)
#' @param grid_n integer grid size if cutpointr not used/available (default 200)
#' @param seed integer or NULL passed to compute_arc_threshold() EM initialization
#' @param verbose logical; print progress messages
#'
#' @return list with:
#'   \item{data}{filtered (or original) df}
#'   \item{applied_mode}{the applied mode ("derived_value","user_value","none")}
#'   \item{arc_info}{list from compute_arc_threshold() (only for derived_value), else NULL}
#'   \item{value_info}{list from compute_value_threshold_via_f1() for derived_value; a stub for user_value; NULL for none}
#'   \item{n_kept}{number of unique IDs retained}
#'   \item{n_total}{number of unique IDs in original df}
#'   \item{thresholds}{list with value_threshold actually used; also arc_threshold for derived_value}
#' @export
apply_reliability_filter <- function(
    df,
    individual_coeffs,
    id_name = "ID",
    time_name = "Time",
    value_name = "Marker_level",
    slope_col = "Slope",
    mode = c("derived_value", "user_value", "none"),
    value_threshold = NULL,
    lower_comp_quantile = 0.99,
    include_negative = TRUE,
    use_cutpointr = TRUE,
    grid_n = 200L,
    seed = NULL,
    verbose = FALSE
) {
  mode <- match.arg(mode)
  
  stopifnot(is.data.frame(df), is.data.frame(individual_coeffs))
  if (!all(c(id_name, slope_col) %in% names(individual_coeffs))) {
    stop("`individual_coeffs` must contain columns: ", id_name, " and ", slope_col)
  }
  if (!all(c(id_name, time_name, value_name) %in% names(df))) {
    stop("`df` must contain columns: ", paste(c(id_name, time_name, value_name), collapse = ", "))
  }
  
  n_total <- length(unique(df[[id_name]]))
  arc_info <- NULL
  value_info <- NULL
  thresholds <- list(value_threshold = NULL, arc_threshold = NULL)
  
  filter_by_baseline_cut <- function(data_long, cut) {
    df_local <- data_long[, c(id_name, time_name, value_name)]
    oo <- order(df_local[[id_name]], df_local[[time_name]])
    df_local <- df_local[oo, , drop = FALSE]
    first_rows <- !duplicated(df_local[[id_name]])  # earliest by sort
    baseline <- df_local[first_rows, c(id_name, value_name), drop = FALSE]
    names(baseline)[2] <- "baseline_value"
    
    keep_ids <- baseline[[id_name]][baseline$baseline_value >= cut]
    data_long[data_long[[id_name]] %in% keep_ids, , drop = FALSE]
  }
  
  if (mode == "none") {
    if (verbose) message("No reliability filter applied (mode='none').")
    return(list(
      data = df,
      applied_mode = mode,
      arc_info = NULL,
      value_info = NULL,
      n_kept = n_total,
      n_total = n_total,
      thresholds = thresholds
    ))
  }
  
  if (mode == "user_value") {
    if (is.null(value_threshold) || !is.numeric(value_threshold) || length(value_threshold) != 1L) {
      stop("For mode='user_value', please supply a single numeric `value_threshold` (same units as ", value_name, ").")
    }
    df_kept <- filter_by_baseline_cut(df, as.numeric(value_threshold))
    n_kept <- length(unique(df_kept[[id_name]]))
    thresholds$value_threshold <- as.numeric(value_threshold)
    if (verbose) message(sprintf("Baseline-value filter (user cut=%.6f): kept %d/%d IDs (%.1f%%).",
                                 thresholds$value_threshold, n_kept, n_total, 100 * n_kept / n_total))
    return(list(
      data = df_kept,
      applied_mode = mode,
      arc_info = NULL,
      value_info = list(
        value_threshold = thresholds$value_threshold,
        summary = NULL,
        baseline_table = NULL,
        warnings = character(0)
      ),
      n_kept = n_kept,
      n_total = n_total,
      thresholds = thresholds
    ))
  }
  
  # mode == "derived_value" ----
  # 1) ARC labels via compute_arc_threshold()
  arc_info <- compute_arc_threshold(
    individual_coeffs = individual_coeffs,
    id_name = id_name,
    slope_col = slope_col,
    lower_comp_quantile = lower_comp_quantile,
    include_negative = include_negative,
    seed = seed,
    verbose = verbose
  )
  thresholds$arc_threshold <- arc_info$arc_threshold
  
  # 2) Learn baseline-value cutpoint by maximizing F1 vs ARC labels
  value_info <- compute_value_threshold_via_f1(
    df = df,
    arc_labels = arc_info$labels,
    id_name = id_name,
    time_name = time_name,
    value_name = value_name,
    engine = if (isTRUE(use_cutpointr)) "cutpointr" else "grid",
    fallback_to_grid = TRUE,
    grid_n = grid_n,
    verbose = verbose
  )
  
  thresholds$value_threshold <- value_info$value_threshold
  
  # 3) Apply the baseline-value cut at earliest visit
  df_kept <- filter_by_baseline_cut(df, thresholds$value_threshold)
  n_kept <- length(unique(df_kept[[id_name]]))
  if (verbose) {
    message(sprintf("Baseline-value filter (derived cut=%.6f): kept %d/%d IDs (%.1f%%).",
                    thresholds$value_threshold, n_kept, n_total, 100 * n_kept / n_total))
  }
  
  list(
    data = df_kept,
    applied_mode = mode,
    arc_info = arc_info,
    value_info = value_info,
    n_kept = n_kept,
    n_total = n_total,
    thresholds = thresholds
  )
}

