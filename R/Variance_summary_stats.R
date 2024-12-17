#' Variance summary statistics
#'
#' @description
#' Estimates the variance effect of several continuous variables jointly
#'
#' @param qt A numeric vector.
#' @param x A data frame, each column represents a covariate that should be numeric.
#' @param iter_num An integer. Represents the number of iterations performed in the Gauss-Newton algorithm
#' @param eps_param A number. The Gauss-Newton algorithm terminates if the incriment change of all variance estimates is smaller than this number.
#'
#' @returns
#' A list with the following objects:
#' * summary, a dataframe with a variance effect estimate for each variable and summary statistics
#' * chi2, the chi2 statistic obtained by considering all parameteres jointly
#' * df, degrees of freedom for the chi2 statistic
#' * pval, p-value of the model
#' * adjusted_values, a vector with qt values that have been adjusted for both mean and variance effects
#'
#' @examples
#' n_val <- 50000
#' x <- as.data.frame(matrix(0,nrow = n_val, ncol = 4))
#' colnames(x) <- c('A','B','C','D')
#' for(i in 1:4) {
#' x[, i] <- rnorm(n_val)
#' }
#' var_vec <- exp(0.2 * x[, 1] - 0.3 * x[, 4])
#' qt_vec <- rnorm(n_val, 0, sqrt(var_vec))
#' res <- var.summary(qt_vec, x)
#' @export var.summary
var.summary <- function(qt, x, iter_num = 50, eps_param = 1e-10) {
  x_matrix <- as.matrix(x)
  n_subjects <- length(qt)
  n_covars <- ncol(x_matrix)

  #First we mean adjust the trait
  A <- as.data.frame(cbind(qt, x))
  colnames(A)[1] <- 'qt'
  qt_mean_adj <- stats::lm(qt ~ ., data = A)$residuals

  #All parameters model
  a_all <- alpha.multi.est(qt_mean_adj, x_matrix, iter_num = iter_num, eps_param = eps_param)
  v_all <- exp(- (x_matrix %*% a_all))
  SigmaSq_all <- sum(qt_mean_adj^2 * v_all) / n_subjects
  l_all <- - n_subjects * log(SigmaSq_all)
  for(k in 1:n_covars){
    l_all <- l_all - sum(a_all[k] * x_matrix[,k])
  }
  adjustment_scalar <- 1 / exp((x_matrix %*% a_all) / 2)
  qt_var_adj <- qt_mean_adj * adjustment_scalar

  #Null model
  a_null <- rep(0, n_covars)
  v_null <- rep(1, n_subjects)
  SigmaSq_null <- sum(qt_mean_adj^2 * v_null) / n_subjects
  l_null <- - n_subjects * log(SigmaSq_null)

  #Comparing all with null
  X2_all_vs_null <- l_all - l_null
  p_all_vs_null <- stats::pchisq(X2_all_vs_null, n_covars, lower.tail = FALSE)

  #Partial models
  M <- as.data.frame(matrix(0, nrow = n_covars, ncol = 5))
  colnames(M) <- c('Variable','var_effect','se','X2','pval')
  M$Variable <- colnames(x)
  M$var_effect <- a_all
  if(n_covars > 1) {
  for(i in 1:n_covars) {
      x_partial <- x_matrix[,-i]
      initial_guess <- a_all[-i]
      a_partial <- alpha.multi.est(qt_mean_adj, x_partial, iter_num = iter_num, eps_param = eps_param, initial_guess = initial_guess)
      v_partial <- exp(- (x_partial %*% a_partial))
      SigmaSq_partial <- sum(qt_mean_adj^2 * v_partial) / n_subjects
      l_partial <- - n_subjects * log(SigmaSq_partial)
      for(k in 1:(n_covars-1)){
        l_partial <- l_partial - sum(a_partial[k] * x_partial[,k])
      }
      X2_partial <- l_all - l_partial
      p_partial <- stats::pchisq(X2_partial, 1, lower.tail = FALSE)
      se_partial <- abs(a_all[i])/sqrt(X2_partial)
      M$se[i] <- se_partial
      M$X2[i] <- X2_partial
      M$pval[i] <- p_partial
    }
  } else {
    M$X2[1] <- X2_all_vs_null
    M$pval[1] <- stats::pchisq(X2_all_vs_null, 1, lower.tail = FALSE)
    M$se[1] <- abs(a_all)/sqrt(X2_all_vs_null)
  }
  return(list(summary = M, chi2 = X2_all_vs_null, df = n_covars, pval = p_all_vs_null, adjusted_values = qt_var_adj))
}
