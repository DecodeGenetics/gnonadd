#' Variant-Environmental interaction effects on a case control variable
#'
#' @description
#' This function estimates the interaction effect of a genetic variant with an environmental factor on a case control variable
#' We apply a simple logistic regression model to estimate interaction effects.
#' We include a linear term for the variant and environmental variable seperately.
#' The variant is coded as 0,1 and 2 for non-carriers, heterozygotes and homozygous carriers of the effect allele.
#' The environmental variable is rank-normalized automatically as part of the function.
#' The interaction term is defined as the product of the genetic and the (normalized) environmental variable.
#' Effect size and significance is based on the interaction term.
#'
#' @param cc A numeric vector
#' @param g A vector with (possibly imputed) genotype values. All entries should be larger than 0 and smaller than 2.
#' @param env A numeric vector with an environmental variable
#' @param yob A numerical vector containing year of birth. If some are unknown they should be marked as -1
#' @param sex A numerical vector containing sex, coded 0 for males, 1 for females and -1 for unknown
#' @param round_imputed A boolian variable determining whether imputed genotype values should be rounded to the nearest integer in the analysis.
#' @param dominance_term A boolian variable determining whether a dominance term for the variant should be included as a covariates in the analysis
#' @param square_env A boolian variable determining whether the square of the environmental trait should be included as a covariate in the analysis
#' @param covariates A dataframe containing any other covariates that should be used; one column per covariate
#'
#' @returns
#' A list with the environmental interaction effect and corresponding standard error, t statistic and p-value
#' @examples
#' g_vec <- rbinom(100000, 2, 0.9)
#' env_vec <- round(runif(100000, min = 0, max = 6))
#' cc_vec <- rbinom(100000,1,0.1 * (1.05^g_vec) *
#'           (1.1^env_vec)* (1.1 ^ (g_vec * env_vec)))
#' res <- env_interaction_CC.calc(cc_vec, g_vec, env_vec)
#' @export
env_interaction_CC.calc <- function(cc, g, env,  yob = rep(-1,length(cc)), sex = rep(-1,length(cc)),
                                    round_imputed = FALSE, dominance_term = FALSE, square_env = FALSE, covariates = as.data.frame(matrix(0, nrow = 0, ncol = 0))){
  r <- rank(env)
  env_normal <- stats::qnorm(r / (length(r) + 1))
  if(round_imputed == TRUE) {
    g <- round(g)
  }
  int <- g * env_normal
  if(stats::sd(int) == 0) {
    warning("Interaction undefined. All interaction values are the same.")
    gamma <- NA
    se <- NA
    z <- NA
    p <- NA
  } else {
    no_date <- yob < 0
    yob[!no_date] <- yob[!no_date] - mean(yob[!no_date])
    yob[no_date] <- 0
    sex <- as.factor(sex)

    #We define a dataframe containing all variables that should be considered
    Env_int_data <- as.data.frame(cbind(cc, int))
    Env_int_data <- cbind(Env_int_data, g)
    Env_int_data <- cbind(Env_int_data, env_normal)
    if(dominance_term == TRUE){
      Env_int_data$g_dom <- as.numeric(round(Env_int_data$g) == 2)
    }
    if(square_env == TRUE){
      Env_int_data$env_square <- (Env_int_data$env_normal)^2
    }
    if(stats::sd(yob) > 0) {
      Env_int_data <- cbind(Env_int_data, yob)
    }
    if(length(unique(no_date)) > 1) {
      Env_int_data <- cbind(Env_int_data, no_date)
    }
    if(length(unique(sex)) > 1) {
      Env_int_data <- cbind(Env_int_data, sex)
    }
    if(nrow(covariates) > 0) {
      Env_int_data <- cbind(Env_int_data, covariates)
    }

    #We use logistic regression to estimate the environmental interaction effect
    l_interaction <- stats::glm(cc ~ ., data = Env_int_data, family = 'binomial')
    param <- "int"
    if(param %in% rownames(stats::coef(summary(l_interaction)))){
      gamma <- summary(l_interaction)$coeff[param, 1]
      se <- summary(l_interaction)$coeff[param, 2]
      z <- summary(l_interaction)$coeff[param, 3]
      p <- summary(l_interaction)$coeff[param, 4]
    }else{
      warning("Singular model matrix")
      gamma <- NA
      se <- NA
      z <- NA
      p <- NA
    }
  }
  return(list(interaction_effect = gamma, standard_error = se, z = z, pval = p))
}
