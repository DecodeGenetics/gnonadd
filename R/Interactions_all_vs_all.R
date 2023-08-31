#' Pairwise interaction effects
#'
#' @description
#' Given a set of variants and a quantitative trait, this function calculates the interaction effect of all possible variant-variant pairs
#'
#'
#' @param qt A numeric vector
#' @param g A matrix, where each colomn represents a variant
#' @param round_imputed A boolian variable determining whether imputed genotype values should be rounded to the nearest integer in the analysis.
#' @param dominance_terms A boolian variable determining whether dominance terms for the variants should be included as covariates in the analysis
#' @param covariates A dataframe containing any other covariates that should be used; one column per covariate
#' @param variant_names A list of the names of the variants
#'
#' @returns
#' A dataframe with all possible variant pairs and their estimated interaction effect
#' @examples
#' g_vec <- matrix(0, nrow = 100000, ncol = 5)
#' freqs <- runif(ncol(g_vec), min = 0, max = 1)
#' for(i in 1:ncol(g_vec)){
#'  g_vec[,i] <- rbinom(100000, 2, freqs[i])
#' }
#'
#' qt_vec <- rnorm(100000) + 0.1 * g_vec[, 1] + 0.2 *
#'           g_vec[, 2] -0.1 * g_vec[, 3] + 0.2 *
#'           g_vec[, 1] * g_vec[, 2]
#' res <- pairwise_int.calc(qt_vec, g_vec)
#' @export
pairwise_int.calc <- function(qt, g, round_imputed = FALSE, dominance_terms = FALSE,
                              variant_names = paste(rep('variant', ncol(g)), as.character(1:ncol(g)), sep="_"),
                              covariates = as.data.frame(matrix(0, nrow = 0, ncol = 0))){
  variant_number <- ncol(g)
  pair_number <- variant_number*(variant_number-1)/2
  A <- data.frame(matrix(0,nrow=pair_number,ncol=6))
  colnames(A) <- c('name1','name2','int_effect','se','t','pval')
  counter <- 0
  for(i in 1 : (variant_number - 1)) {
    for(j in (i+1):variant_number){
      counter <- counter + 1
      A$name1[counter] <- variant_names[i]
      A$name2[counter] <- variant_names[j]
      res <- interaction.calc(qt, g[, i], g[, j], round_imputed = round_imputed,
                              dominance_terms = dominance_terms, covariates = covariates)
      A$int_effect[counter] <- res$interaction_effect
      A$se[counter] <- res$standard_error
      A$t[counter] <- res$t
      A$pval[counter] <- res$pval

    }
  }
  return(A)
}
