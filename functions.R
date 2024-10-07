#' Generates synthetic data
#'
#'
#' @param n Number of individuals
#' @param T total length
#' @param l length of time between observations
#' @return all_data synthetic dataset with n participants with unique ID's, observations at an interval of Time l, and Marker_level values based on distribution parameters
#' @export

# Generate subsets of individual disease progression data using 5PL
generate_synth_data <- function (n, T, l){
  b<- runif(n, 10,16)
  a <- rtruncnorm(n, 0, Inf, 0.4, 0.2)
  c <- rtruncnorm(n, 0, Inf, 0.01, 0.005)
  f <- rtruncnorm(n, -Inf, 1, 0.98, 0.005)
  g <- rnorm(n, 20, 1)
  b_star <- b+log(2^(1/g)-1)/a
  tau_all <- seq(0, T, by = l)
  
  all_data <- data.frame()
  # generate data for each individual 
  for (i in 1: n){
    R_tau <- c[i] + (f[i]-c[i])/(1+exp(-a[i]*(tau_all-b_star[i])))^g[i]
    # combing individual data into the data frame
    individual_data<- data.frame(ID = i, Time = tau_all, Marker_level = R_tau)
    all_data <- rbind(all_data, individual_data)
  }
  return(all_data)
}