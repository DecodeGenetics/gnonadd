#' Pairwise environmental interaction effects
#'
#' @description
#' Given a set of variants and environmental traits, and a single quantitative trait, this function calculates the interaction effect of all possible variant-environmental pairs
#'
#'
#' @param qt A numeric vector
#' @param g A matrix, where each colomn represents a variant
#' @param env A matrix, where each row represents an environmental variable
#' @param round_imputed A boolian variable determining whether imputed genotype values should be rounded to the nearest integer in the analysis.
#' @param dominance_term A boolian variable determining whether a dominance term for the variant should be included as a covariates in the analysis
#' @param square_env A boolian variable determining whether the square of the environmental trait should be included as a covariate in the analysis
#' @param covariates A dataframe containing any other covariates that should be used; one column per covariate
#' @param variant_names A list of the names of the variants
#' @param env_names A list of the names of the environmental variables
#'
#' @returns
#' A dataframe with all possible variant-environmental pairs and their estimated interaction effect
#' @examples
#' g_vec <- matrix(0, nrow = 100000, ncol = 3)
#' freqs <- runif(ncol(g_vec), min = 0, max = 1)
#' env_vec <- matrix(0, nrow = 100000, ncol = 3)
#' for(i in 1:ncol(g_vec)){
#'  g_vec[, i] <- rbinom(100000, 2, freqs[i])
#' }
#' for( i in 1:ncol(env_vec)){
#'  env_vec[, i] <- round(runif(100000,min=0,max=6))
#' }
#'
#' qt_vec <- rnorm(100000) + 0.1 * g_vec[, 1] + 0.2 *
#'           g_vec[, 2] -0.1 * env_vec[, 3] + 0.1 *
#'           env_vec[, 1] + 0.1 * g_vec[, 1] * env_vec[, 1]
#' res <- pairwise_env_int.calc(qt_vec, g_vec, env_vec)
#' @export
pairwise_env_int.calc <- function(qt, g, env, round_imputed = FALSE, dominance_term = FALSE, square_env = FALSE,
                                  covariates = as.data.frame(matrix(0, nrow = 0, ncol = 0)),
                                  variant_names = paste(rep('variant', ncol(g)), as.character(1:ncol(g)), sep="_"),
                                  env_names =  paste(rep('env', ncol(env)), as.character(1:ncol(env)), sep="_")){

  pair_number <- ncol(g) * ncol(env)
  A <- data.frame(matrix(0, nrow = pair_number, ncol=6))
  colnames(A) <- c('variant_name', 'env_name', 'int_effect', 'se', 't', 'pval')
  counter <- 0
  for(i in 1:ncol(g)) {
    for(j in 1:ncol(env)){
      counter <- counter + 1
      A$variant_name[counter] <- variant_names[i]
      A$env_name[counter] <- env_names[j]
      res <- env_interaction.calc(qt, g[, i], env[, j], round_imputed = round_imputed, dominance_term = dominance_term,
                                  square_env = square_env, covariates = covariates )
      A$int_effect[counter] <- res$interaction_effect
      A$se[counter] <- res$standard_error
      A$t[counter] <- res$t
      A$pval[counter] <- res$pval

    }
  }
  return(A)
}



