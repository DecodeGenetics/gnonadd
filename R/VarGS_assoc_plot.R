#' Actual variance vs predicted variance plot
#'
#' @description
#' This tool creates a line plot that compares the predicted variance of data to its actual variance.
#'
#' @param qt A numeric vector.
#' @param v_score A numeric vector.
#' @param bins An integer.
#' @param xlab A string.
#' @param ylab A string.
#' @param title A string.
#'
#' @returns
#' A plot comparing predicted variance to actual variance.
#'
#' @examples
#' n_val <- 100000L
#' v_vec <- exp(rnorm(n_val, 0, 0.1))
#' qt_vec <- stats::rnorm(n_val, 0, sqrt(v_vec))
#' VarGS.plot(qt_vec, v_vec)
#' @export
VarGS.plot <- function(qt, v_score, bins = 10, xlab = 'Predicted variance', ylab = 'Variance', title = '') {
  vRS_mean <- NULL
  quantile_variance <- NULL
  n <- length(qt)
  v_score <- log(v_score)
  v_score_quantiles <- as.numeric(stats::quantile(v_score, 0 : bins / bins))
  quantile_list <- rep(0,n)
  for(i in 1:bins){
    quantile_list[v_score > v_score_quantiles[i] & v_score <= v_score_quantiles[i + 1]] <- i
  }
  M <- as.data.frame(matrix(0, bins, 3))
  colnames(M) <- c('quantile', 'vRS_mean', 'quantile_variance')
  M$quantile <- 1:bins
  for(i in 1:bins){
    M$vRS_mean[i] <- mean(v_score[quantile_list == i])
    M$quantile_variance[i] <- stats::var(qt[quantile_list == i])
  }
  M$vRS_mean <- exp(M$vRS_mean)
  ggplot2::ggplot(M, ggplot2::aes(x = vRS_mean, y = quantile_variance)) +
    ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::theme_classic() +
    ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + ggplot2::ggtitle(title)
}
