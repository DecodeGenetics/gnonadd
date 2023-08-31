#' Genetic dominance effects on a case control variable
#'
#' @description
#' This function estimates the dominance effect of a genetic variant on a case-control variable
#' We apply a logistic regression model to estimate dominance effects.
#' We include a linear term, coded as 0,1 and 2 for non-carriers, heterozygotes and homozygous carriers of the effect allele.
#' We also include a dominance term, coded as 1 for homozygous carriers and 0 for others.
#' Effect size and significance is based on the dominance term.
#'
#' @param cc A case control vector, containing 0's and 1's
#' @param g A vector with (possibly imputed) genotype values. All entries should be larger than 0 and smaller than 2.
#' @param yob A numerical vector containing year of birth. If some are unknown they should be marked as -1
#' @param sex A numerical vector containing sex, coded 0 for males, 1 for females and -1 for unknown
#' @param round_imputed A boolian variable determining whether imputed genotype values should be rounded to the nearest integer in the analysis
#' @param covariates A dataframe containing any other covariates that should be used; one column per covariate.
#'
#' @returns
#' A list with the dominanc effect (on log-scale) and corresponding standard error, z statistic and p-value
#' @examples
#' g_vec <- rbinom(100000, 2, 0.3)
#' cc_vec <- rbinom(100000, 1, 0.1 * (1.2 ^ (g_vec^2)))
#' res <- dominance_CC.calc(cc_vec, g_vec)
#' @export
dominance_CC.calc <- function(cc, g, yob=rep(-1,length(cc)), sex=rep(-1,length(cc)), round_imputed = FALSE, covariates = as.data.frame(matrix(0, nrow = 0, ncol = 0))){
  g_rounded <- round(g)
  if(round_imputed == TRUE){
    g <- round(g)
  }
  if(length(unique(as.factor(g_rounded))) < 3) {
    warning("Dominance effect undefined. There are no subjects of one or more genotype group.")
    delta <- NA
    se <- NA
    z <- NA
    p <- NA
  } else {
    no_date <- yob < 0
    yob[!no_date] <- yob[!no_date] - mean(yob[!no_date])
    yob[no_date] <- 0
    sex <- as.factor(sex)
    g <- g - mean(g)
    g2 <- as.numeric(g_rounded == 2)

    #We define a dataframe containing all variables that should be considered
    Dom_data <- as.data.frame(cbind(cc, g2))
    Dom_data <- cbind(Dom_data, g)
    if(stats::sd(yob) > 0) {
      Dom_data <- cbind(Dom_data, yob)
    }
    if(length(unique(no_date)) > 1) {
      Dom_data <- cbind(Dom_data, no_date)
    }
    if(length(unique(sex)) > 1) {
      Dom_data <- cbind(Dom_data, sex)
    }
    if(nrow(covariates) > 0) {
      Dom_data <- cbind(Dom_data, covariates)
    }

    #We use logistic regression to estimate the dominance effect
    l_delta <- stats::glm(cc ~ ., data = Dom_data, family = 'binomial')
    param <- "g2"
    if(param %in% rownames(stats::coef(summary(l_delta)))){
      delta <- summary(l_delta)$coeff[param, 1]
      se <- summary(l_delta)$coeff[param, 2]
      z <- summary(l_delta)$coeff[param, 3]
      p <- summary(l_delta)$coeff[param, 4]
    }else{
      warning("Singular model matrix")
      delta <- NA
      se <- NA
      z <- NA
      p <- NA
    }
  }
  return(list(dominance_effect = delta, standard_error = se, z = z, pval = p))
}

