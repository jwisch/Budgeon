% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrapping_function.R
\name{bootstrap_get_Time_to_Positivity}
\alias{bootstrap_get_Time_to_Positivity}
\title{Wrapper function for bootrapping the model fits}
\usage{
bootstrap_get_Time_to_Positivity(
  df,
  PET_pos_threshold,
  id_name,
  time_name,
  value_name,
  num_bootstraps = 1000,
  bootstrap_percent = 0.8,
  degree = 3,
  printIter = TRUE
)
}
\arguments{
\item{df}{dataframe}

\item{PET_pos_threshold}{numeric threshold for biomarker positivity}

\item{id_name}{column name that contains the subject id value, e.g. "ID"}

\item{time_name}{column name that contains the time value, e.g. "Time"}

\item{value_name}{column name that contains the biomarker value, e.g. "Value"}

\item{bootstrap_percent}{percentage of IDs to include in each bootstrap iteration}

\item{degree}{degree of polynomial fit}

\item{printIter}{print the iterations of the bootstrap as a mode of progress}

\item{num_bootraps}{number of bootstraps}
}
\value{
df_tau dataframe containing ID, predicted value based on model fit, time to positivity estimate
}
\description{
Wrapper function for bootrapping the model fits
}
