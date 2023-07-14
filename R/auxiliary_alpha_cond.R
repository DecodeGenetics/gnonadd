#' Two auxiliary functions for the alpha.cond function
#'
#' @description
#' G_likely is a function with the same zeros as the gradiant of the likelyhood function corresponding to the alpha.cond analysis.
#' Hess_likely is the Jacobian of G_likely, thus it serves the same purpose as the Hessian of the likelyhood function corresponding to the alpha.cond analysis.
#'
#' @param a A numeric vector consisting of variance effect estimates for each variant being tested.
#' @param rss An array containing the residual sum of squares for each possible genotype combination for the variants being tested.
#' @param freqs A numeric vector containing the frequencies of the variants being tested.
#' @noRd
G_likely <- function(a, rss, freqs) {
  n <- length(a)
  vecs <- list(c(0, 1, 2))   #possible genotypes for a single variant
  vecs <- rep(vecs, n)      #possible genotypes for all variants
  vecs <- as.matrix(do.call(expand.grid, vecs))   #All possible genotype combinations
  dl_da <- rep(0, n)
  for(j in 1:n){
    for(i in 1:nrow(vecs)){
      dl_da[j] <- dl_da[j]+(2 * freqs[j] - vecs[i, j]) * rss[t(1 + vecs[i, ])] * prod(a^(2 - vecs[i, ]))
    }
  }
  return(dl_da)
}
Hess_likely <- function(a, rss, freqs) {
  n <- length(a)
  vecs <- list(c(0, 1, 2))   #possible genotypes for a single variant
  vecs <- rep(vecs, n)      #possible genotypes for all variants
  vecs <- as.matrix(do.call(expand.grid, vecs))   #All possible genotype combinations
  H <- matrix(0, n, n)
  for(r in 1:n){
    for(s in 1:n){
      temp_s <- 1:n
      unit_s <- as.numeric(temp_s == s)
      for(i in 1:nrow(vecs)){
        H[r, s] <- H[r, s]+(2 * freqs[r]-vecs[i, r])*(2 - vecs[i, s]) * rss[t(1+vecs[i, ])] * prod(a^(2 - vecs[i,] - unit_s))
      }
    }
  }
  return(H)
}
