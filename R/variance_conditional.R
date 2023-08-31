#' Conditional analysis for genetic variance effects
#'
#'
#' @description
#' We estimate the variance effect of a primary variant conditioned on one or more secondary variants.
#' We apply a likelyhood ratio test with one degree of freedom
#'
#' H0: All secondary variants have a variance effect, but not the primary one.
#' H1: All variants have a variance effect, including the primary one
#'
#' Under both models we assume a full mean effect model. That is, the number of mean-value parameters is 3^(n+1), where n is the number of covariates
#' Thus the null model has 3^(n+1)+n degrees of freedom whereas the alternative model has 3^(n+1)+n+1
#' Due to the exponential growth of free parameters, the test might have low statistical power if many covariates are used
#'
#' IMPORTANT NOTE: We use the Gauss-Newton algorithm to estimate many parameters.
#' We do not check if the algorithm converges to the true minimum values, or if it converges at all.
#' Thus, we advise against blindly believing all results
#'
#' @param qt A numeric vector.
#' @param g An integer vector.
#' @param g_covar An integer matrix where each column corresponds to a genetic covariate
#' @param iter_num An integer. Represents the number of iterations performed in the Gauss-Newton algorithm
#'
#' @returns
#' A list with the values:
#' * alpha, the estimated variance effect, conditioned on the covariates
#' * pval, the p-value corresponding to alpha

#' @examples
#' n <- 10000
#' qt <- rnorm(n)
#' g <- rbinom(n, 2, 0.3)
#' qt <- qt * 1.2^g
#' g_covar <- as.data.frame(matrix(0, nrow = n, ncol = 3))
#' for(i in 1:ncol(g_covar)){
#'   freq <- runif(1, min = 0, max = 1)
#'   g_covar[, i] <- rbinom(n, 2, freq)
#' }
#' res <- alpha.cond(qt, g, g_covar)
#' @export
alpha.cond <- function(qt, g, g_covar, iter_num=50) {
  k <- length(qt)    #Number of measurements
  n <- ncol(g_covar) + 1  #Number of variants (counting the main variant)
  A <- cbind(g, g_covar)
  flips <- rep(1, ncol(A)) #We flip to the minor allele of each variant for techincal reasons. We keep track of which ones were flipped and flip back in the end.
  for(i in 1:ncol(A)){
    A[, i] <- round(A[, i])
    if(mean(A[, i]) > 1){
      flips[i] <- -1
      A[, i] <- 2 - A[, i]
    }

  }
  freqs <- rep(0, n) #A vector with the frequencies of the variants
  for(i in 1:ncol(A)){
    freqs[i] <-  sum(A[, i]) / (2 * k)
  }
  vecs <- list(c(0, 1, 2))   #possible genotypes for a single variant
  vecs <- rep(vecs, n)      #possible genotypes for all variants
  vecs <- as.matrix(do.call(expand.grid, vecs))   #All possible genotype combinations

  # The followign if-loop is here to save running time. The fastest code depends on how many covariates there are.
  # If there are only a few covariates, it makes more sense to loop over all possible genotype groops
  # If there are more covariates, it makes more sense to loop over all the subjects and check what genotype groop each one belongs to
  # We tested the code with n=300.000 and with that datasize, the cutoff seems to be around 5.
  covar_cutoff <- 5
  if(n > covar_cutoff){
    m <- array(0, rep(3, n))  #An array storing the mean value for each possible genotype-group
    number_counter <- array(0, rep(3, n)) #An array storing the number of subjects in each genotype-group
    rss <- array(0, rep(3, n)) #An array storing the residual sum of squares
    A_temp <- as.matrix(A)
    for(i in 1:k){
      m[t(1 + A_temp[i, ])] <- m[t(1 + A_temp[i, ])]+qt[i]
      number_counter[t(1 + A_temp[i, ])] <- number_counter[t(1 + A_temp[i,])] + 1
    }
    m <- m / number_counter
    for(i in 1:k){
      rss[t(1 + A_temp[i, ])] <- rss[t(1 + A_temp[i, ])] + (qt[i] - m[t(1 + A_temp[i, ])])^2
    }
  }else{
    m <- array(0, rep(3, n)) #An array storing the mean value for each possible genotype-group
    rss <- array(0, rep(3, n)) #An array storing the residual sum of squares
    for(i in 1:nrow(vecs)){
      subject_subset <- rowSums(t(vecs[i, ] == t((A[, 1:ncol(A)])))) == n  #The list of subjects with genotype group equal to vecs[i,]
      if(sum(subject_subset) == 0){
        m[t(1 + vecs[i,])] <- 0
        rss[t(1 + vecs[i,])] <- 0
      }else{
        m[t(1 + vecs[i,])] <- mean(qt[subject_subset])
        rss[t(1 + vecs[i,])] <- sum((qt[subject_subset] - m[t(1 + vecs[i, ])])^2)
      }
    }
  }

  #Null model.
  #The problem can be solved analytically if there is only one covariate. Hence the if-statement
  if(n <= 2){
    rss_null <- c(0, 0, 0)
    rss_null[1] <- rss[1, 1]+rss[2, 1]+rss[3, 1]
    rss_null[2] <- rss[1, 2]+rss[2, 2]+rss[3, 2]
    rss_null[3] <- rss[1, 3]+rss[2, 3]+rss[3, 3]
    aa <- 2 * freqs[2] * rss_null[1]
    bb <- (2 * freqs[2] - 1) * rss_null[2]
    cc <- (2 * freqs[2] - 2) * rss_null[3]
    alpha_null <- (- bb + sqrt(bb^2 - 4 * aa * cc)) / (2 * aa)    #This is the maximum likelyhood estimator for alpha
    sigma_null <- (rss_null[3] / (alpha_null^2) + rss_null[2] / alpha_null + rss_null[1]) / k  #This is the maximum likelyhood estimator for sigma
    l_null <- - k * log(sigma_null) / 2 - log(alpha_null) * k * freqs[2]   #This is the likelyhood function evaluated at those values
  }else{

    #if there is more than one covariate we use Newtons method to solve the problem numerically
    vecs_null <- list(c(0, 1, 2))
    vecs_null <- rep(vecs_null, n - 1)
    vecs_null <- as.matrix(do.call(expand.grid, vecs_null))
    rss_null <- array(0, rep(3, n - 1))
    for(i in 1:nrow(vecs_null)){
      rss_null[t(1 + vecs_null[i, ])] <- rss[t(1 + c(0, vecs_null[i, ]))] + rss[t(1 + c(1, vecs_null[i, ]))] + rss[t(1 + c(2, vecs_null[i, ]))]
    }
    freqs_null <- freqs[2:n]
    alpha_null <- rep(1, n - 1)  # This is the orignial guess of the alpha parameters. The subsequent for-loop is Newtons method to obtain the correct maximum likelyhood estimators
    for(i in 1:iter_num){
      alpha_null <- tryCatch(alpha_null - solve(Hess_likely(alpha_null, rss_null, freqs_null)) %*% G_likely(alpha_null, rss_null, freqs_null), error = function(err) NA)
    }
    sigma_null <-0
    for(i in 1:nrow(vecs_null)){
      sigma_null <- sigma_null + rss_null[t(1 + vecs_null[i, ])] * prod(alpha_null^(2 - vecs_null[i, ]))
    }
    sigma_null <- sigma_null / (k * prod(alpha_null^2)) # This is the maximum likelyhood estimator for sigma
    l_null <- - k * log(sigma_null) / 2 - k * sum(log(alpha_null) * freqs_null)  ### This is the likelyhood function evaluated at those values
  }

  #Alt model
  #We always need newtons method to solve for the parameters in the alternative model.
  alpha_alt <- rep(1, n)
  for(i in 1:iter_num){
    alpha_alt <- tryCatch(alpha_alt - solve(Hess_likely(alpha_alt, rss, freqs)) %*% G_likely(alpha_alt, rss, freqs), error = function(err) NA)
  }
  if(is.na(alpha_alt[1])){
    P <- alpha_alt
  }else{
    sigma_alt <- 0
    for(i in 1:nrow(vecs)){
      sigma_alt <- sigma_alt + rss[t(1 + vecs[i, ])] * prod(alpha_alt^(2 - vecs[i, ]))
    }
    sigma_alt <- sigma_alt / (k * prod(alpha_alt^2)) # This is the maximum likelyhood estimator for sigma
    l_alt <- - k * log(sigma_alt) / 2 - k * sum(log(alpha_alt)*freqs)  ### This is the likelyhood function evaluated at those values

    #Finally we compare the maximized likelyhood functions to obtain an approximated X2-statistic
    X2 <- 2 * (l_alt - l_null)
    alpha_alt <- alpha_alt^(flips)
    P <- stats::pchisq(X2, 1, lower.tail = FALSE)
  }
  return(list(alpha = alpha_alt[1], pval = P))
}
