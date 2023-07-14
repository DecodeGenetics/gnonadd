#' Expected variance effect from additive effect
#'
#' @description
#' This function interpolates data from a simple simulation to give an estimate
#' of the variance effect induced by an additive effect. The simulation code
#' is stored under inst/raw/. We assume that the trait has been inverse normal
#' transformed. Under the simulation, there is no variance effect, so the variance
#' effect is fully induced by the inverse normal transform.
#'
#' @param maf Minor allele frequency of the variant, should be in the range 0 to 0.5.
#' @param beta_add Additive effect of the variant, should be in the range 0 to 3.5.
#'                 This variable can be a vector of values.
#'
#' @returns
#' The expected variance effect for the variant from the given maf, beta combination.
#' @examples
#' maf <- 0.1
#' beta_val <- 0.3
#' expected_var <- expected.variance.effect(maf, beta_val)
#'
#' beta_vec <- seq(0.1,2, length.out = 20)
#' expected_var <- expected.variance.effect(maf, beta_vec)
#' @export
expected.variance.effect <- function(maf, beta_add) {
  beta_add <- abs(beta_add)
  asymptote_warning <- "Input beta out of range, higher than inverse normal asymptote\nReturning NA for some beta input\n"
  warning_has_appeared <- FALSE
  f_alleles <- sort(unique(.variance.simulation$f_par))
  # Minor input checks
  if(maf < min(f_alleles)) {
    stop(paste0("Cannot currently extrapolate for maf < ", min(f_alleles)))
  }
  if(maf > 0.5) {
    stop("maf should be in the range ]0, 0.5]")
  }

  # If we are this close to an existing maf, we just use it
  eps_param <- min(f_alleles) / 2.0
  beta_len <- length(beta_add)
  alpha_out <- rep(NA, beta_len)
  # Case 1: Close to existing curve
  if(any(abs(f_alleles - maf) < eps_param)) {
    f_use <- which(abs(f_alleles - maf) < eps_param)[1] # Use first in case of ties...
    path_data <- .variance.simulation[.variance.simulation$f_par == f_alleles[f_use],]
    alpha_smooth <- stats::smooth(path_data$alpha_hat)
    for(i in seq_len(beta_len)) {
      # Check if input beta is out of range
      if(beta_add[i] > max(path_data$beta_hat)) {
        if(!warning_has_appeared) {
          warning(asymptote_warning)
          warning_has_appeared <- TRUE
        }
        next
      }
      ind_out <- which.min(abs(path_data$beta_hat - beta_add[i]))
      # Let's just take the closest value
      alpha_out[i] <- alpha_smooth[ind_out]
    }
    return(alpha_out)
  }

  # Case 2: In between maf values
  higher_ind <- min(which(f_alleles > maf))
  lower_ind <- higher_ind - 1

  # Create new interpolated smoothed alpha and beta_hat
  len_interval <- f_alleles[higher_ind] - f_alleles[lower_ind]
  t_val <- (maf - f_alleles[lower_ind]) / len_interval

  path_data_lower <- .variance.simulation[.variance.simulation$f_par == f_alleles[lower_ind],]
  path_data_higher <- .variance.simulation[.variance.simulation$f_par == f_alleles[higher_ind],]
  alpha_smooth_lower <- stats::smooth(path_data_lower$alpha_hat)
  alpha_smooth_higher <- stats::smooth(path_data_higher$alpha_hat)

  # Linear interpolation of the curves
  beta_hat_interp <- (t_val * path_data_higher$beta_hat) + (1 - t_val) * path_data_lower$beta_hat
  alpha_smooth_interp <- (t_val * alpha_smooth_higher) + (1 - t_val) * alpha_smooth_lower

  for(i in seq_len(beta_len)) {
    # Check if input beta is out of range
    if(beta_add[i] > max(beta_hat_interp)) {
      if(!warning_has_appeared) {
        warning(asymptote_warning)
        warning_has_appeared <- TRUE
      }
      next
    }
    ind_out <- which.min(abs(beta_hat_interp - beta_add[i]))
    # Let's just take the closest value
    alpha_out[i] <- alpha_smooth_interp[ind_out]
  }
  return(alpha_out)
}
