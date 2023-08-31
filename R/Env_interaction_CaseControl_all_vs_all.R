#' Pairwise environmental interaction effects for a case control variable
#'
#' @description
#' Given a set of variants and environmental traits, and a single case control variable, this function calculates the interaction effect of all possible variant-environmental pairs
#'
#' @param cc A numeric vector
#' @param g A matrix, where each colomn represents a variant
#' @param env A matrix, where each row represents an environmental variable
#' @param yob A numerical vector containing year of birth. If some are unknown they should be marked as -1
#' @param sex A numerical vector containing sex, coded 0 for males, 1 for females and -1 for unknown
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
#' N_run <- 25000
#' g_vec <- matrix(0, nrow = N_run, ncol = 3)
#' freqs <- runif(ncol(g_vec), min = 0, max = 1)
#' env_vec <- matrix(0, nrow = N_run, ncol = 3)
#' for(i in 1:ncol(g_vec)){
#'  g_vec[, i] <- rbinom(N_run, 2, freqs[i])
#' }
#' for( i in 1:ncol(env_vec)){
#'  env_vec[, i] <- round(runif(N_run,min=0,max=6))
#' }
#' cc_vec <- rbinom(N_run,1,0.1 * (1.05 ^ g_vec[, 1]) *
#'           (1.06 ^ env_vec[,1]) * (0.95 ^ g_vec[, 2]) *
#'           (1.1^(g_vec[, 1] * env_vec[, 1])))
#' res <- pairwise_env_int_CC.calc(cc_vec, g_vec, env_vec)
#' @export
pairwise_env_int_CC.calc <- function(cc, g, env, yob = rep(-1,length(cc)), sex = rep(-1,length(cc)),
                                     round_imputed = FALSE, dominance_term = FALSE, square_env = FALSE, covariates = as.data.frame(matrix(0, nrow = 0, ncol = 0)),
                                     variant_names = paste(rep('variant', ncol(g)), as.character(1:ncol(g)), sep="_"),
                                     env_names =  paste(rep('env', ncol(env)), as.character(1:ncol(env)), sep="_")){
  pair_number <- ncol(g) * ncol(env)
  A <- data.frame(matrix(0, nrow = pair_number, ncol=6))
  colnames(A) <- c('variant_name', 'env_name', 'int_effect', 'se', 'z', 'pval')

  counter <- 0
  for(i in 1:ncol(g)) {
    for(j in 1:ncol(env)){
      counter <- counter + 1
      A$variant_name[counter] <- variant_names[i]
      A$env_name[counter] <- env_names[j]
      res <- env_interaction_CC.calc(cc, g[, i], env[, j], yob = yob, sex = sex,
                                     round_imputed = round_imputed, dominance_term = dominance_term, square_env = square_env, covariates = covariates )
      A$int_effect[counter] <- res$interaction_effect
      A$se[counter] <- res$standard_error
      A$z[counter] <- res$z
      A$pval[counter] <- res$pval

    }
  }
  return(A)
}



