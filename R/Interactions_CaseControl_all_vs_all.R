#' Pairwise interaction effects for a case control variable
#'
#' @description
#' Given a set of variants and a case control variable, this function calculates the interaction effect of all possible variant-variant pairs
#'
#' @param cc A numeric vector
#' @param g A matrix, where each colomn represents a variant
#' @param yob A numerical vector containing year of birth. If some are unknown they should be marked as -1
#' @param sex A numerical vector containing sex, coded 0 for males, 1 for females and -1 for unknown
#' @param round_imputed A boolian variable determining whether imputed genotype values should be rounded to the nearest integer in the analysis.
#' @param dominance_terms A boolian variable determining whether dominance terms for the variants should be included as covariates in the analysis
#' @param covariates A dataframe containing any other covariates that should be used; one column per covariate
#' @param variant_names A list of the names of the variants
#'
#' @returns
#' A dataframe with all possible variant pairs and their estimated interaction effect
#' @examples
#' N_run <- 25000
#' g_vec <- matrix(0, nrow = N_run, ncol = 5)
#' freqs <- runif(ncol(g_vec), min = 0,max = 1)
#' for(i in 1:ncol(g_vec)){
#'  g_vec[, i] <- rbinom(N_run, 2, freqs[i])
#' }
#' cc_vec <- rbinom(N_run,1,0.1 * (1.05 ^ g_vec[, 1]) *
#'           (1.06 ^ g_vec[,2]) * (0.95 ^ g_vec[, 3]) *
#'           (1.5^(g_vec[,1] * g_vec[,2])))
#' res <- pairwise_int_CC.calc(cc_vec, g_vec)
#' @export
pairwise_int_CC.calc <- function(cc, g, yob = rep(-1,length(cc)), sex = rep(-1,length(cc)),
                                 round_imputed = FALSE, dominance_terms = FALSE, covariates = as.data.frame(matrix(0, nrow = 0, ncol = 0)),
                                 variant_names = paste(rep('variant', ncol(g)), as.character(1:ncol(g)), sep="_")){
   variant_number <- ncol(g)
  pair_number <- variant_number*(variant_number-1)/2
  A <- data.frame(matrix(0,nrow=pair_number,ncol=6))
  colnames(A) <- c('name1','name2','int_effect','se','z','pval')
  counter <- 0
  for(i in 1 : (variant_number - 1)) {
    for(j in (i+1):variant_number){
      counter <- counter + 1
      A$name1[counter] <- variant_names[i]
      A$name2[counter] <- variant_names[j]
      res <- interaction_CC.calc(cc, g[, i], g[, j], yob = yob, sex = sex, round_imputed = round_imputed,
                                 dominance_terms = dominance_terms, covariates = covariates)
      A$int_effect[counter] <- res$interaction_effect
      A$se[counter] <- res$standard_error
      A$z[counter] <- res$z
      A$pval[counter] <- res$pval

    }
  }
  return(A)
}
