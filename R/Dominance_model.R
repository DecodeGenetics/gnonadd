#' Genetic dominance effects
#'
#' @description
#' This function estimates the dominance effect of a genetic variant on a quantitatvie trait
#' Nothing fancy here. We apply a simple linear regression model to estimate dominance effects.
#' We include a linear term, coded as 0,1 and 2 for non-carriers, heterozygotes and homozygous carriers of the effect allele.
#' We also include a dominance term, coded as 1 for homozygous carriers and 0 for others.
#' Effect size and significance is based on the dominance term.
#'
#' @param qt A numeric vector
#' @param g A vector with (possibly imputed) genotype values. All entries should be larger than 0 and smaller than 2.
#' @param round_imputed A boolian variable determining whether imputed genotype values should be rounded to the nearest integer in the analysis
#' @param covariates A dataframe containing any covariates that should be used; one column per covariate.
#'
#' @returns
#' A list with the dominanc effect and corresponding standard error, t statistic and p-value
#' @examples
#' g_vec <- rbinom(100000, 2, 0.3)
#' qt_vec <- rnorm(100000) + 0.2 * g_vec^2
#' res <- dominance.calc(qt_vec, g_vec)
#' @export
dominance.calc <- function(qt, g, round_imputed = FALSE, covariates = as.data.frame(matrix(0, nrow = 0, ncol = 0))){
  g_rounded <- round(g)
  if(round_imputed == TRUE){
    g <- round(g)
  }
  if(length(unique(as.factor(g_rounded))) < 3) {
    warning("Dominance effect undefined. There are no subjects of one or more genotype group.")
    delta <- NA
    se <- NA
    t <- NA
    p <- NA
  }else {
    g2 <- as.numeric(g_rounded == 2)

    #We define a dataframe containing all variables that should be considered
    Dom_data <- as.data.frame(cbind(qt, g2))
    Dom_data <- cbind(Dom_data, g)
    if(nrow(covariates) > 0) {
      Dom_data <- cbind(Dom_data, covariates)
    }

    #We use linear regression to estimate the dominance effect
    l_delta <- stats::lm(qt ~ ., data = Dom_data)
    param <- "g2"
    if(param %in% rownames(stats::coef(summary(l_delta)))){
      delta <- summary(l_delta)$coeff[param, 1]
      se <- summary(l_delta)$coeff[param, 2]
      t <- summary(l_delta)$coeff[param, 3]
      p <- summary(l_delta)$coeff[param, 4]
    }else{
      warning("Singular model matrix")
      delta <- NA
      se <- NA
      t <- NA
      p <- NA
    }
  }
  return(list(dominance_effect = delta, standard_error = se, t = t, pval = p))
}
