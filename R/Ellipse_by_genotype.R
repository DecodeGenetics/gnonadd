#' Ellipse best fit plot
#'
#' @description
#' This tool creates a scatter plot along with regression lines. Additionally it finds and plots the best ellipses that fit the data.
#'
#' @param qt1 A numeric vector.
#' @param qt2 A numeric vector.
#' @param g An integer vector.
#' @param trait_name1 A string.
#' @param trait_name2 A string.
#' @param title A string.
#' @param sample_size A positive integer.
#' @returns
#' A scatter plot.
#' @examples
#' n_val <- 10000L
#' geno_vec <- c(rep(0, n_val), rep(1, n_val), rep(2, n_val))
#' qt_g0 <- MASS::mvrnorm(n_val, mu = c(0, 0), Sigma = matrix(c(0.93, 0.88, 0.88, 0.92), ncol = 2))
#' qt_g1 <- MASS::mvrnorm(n_val, mu = c(0, 0), Sigma = matrix(c(0.98, 0.88, 0.88, 0.90), ncol = 2))
#' qt_g2 <- MASS::mvrnorm(n_val, mu = c(0, 0), Sigma = matrix(c(1.57, 0.81, 0.81, 0.59), ncol = 2))
#' qt_vec <- rbind(qt_g0, qt_g1)
#' qt_vec <- rbind(qt_vec, qt_g2)
#' res <- ellipse.by.gen(qt_vec[, 1], qt_vec[, 2], geno_vec)
#' @export ellipse.by.gen
ellipse.by.gen <- function(qt1, qt2, g, trait_name1 = 'qt trait 1', trait_name2 = 'qt trait 2',
                        title = '', sample_size = 500) {
  g_factor <- NULL
  x0 <- NULL
  x1 <- NULL
  x2 <- NULL
  y0 <- NULL
  y1 <- NULL
  y2 <- NULL
  g <- round(g)
  D <- cbind(qt1, qt2)
  D <- cbind(D, g)
  D <- as.data.frame(D)
  colnames(D) <- c('qt1', 'qt2', 'g')
  D$g_factor <- factor(D$g, levels = 0:2, labels = c('Non-carriers', 'Heterozygotes', 'Homozygotes'))
  M <- as.data.frame(matrix(0,500,7))
  colnames(M) <- c('t','x0','y0','x1','y1','x2','y2')
  M$t <- (1:500/500)*2*pi
  D_sample <- D[c(),]
  Arrow_data <- as.data.frame(matrix(0,6,4))
  colnames(Arrow_data) <- c('start_x', 'start_y', 'end_x', 'end_y')
  for(i in 0:2) {
    D_temp <- D[D$g == i, ]
    if(nrow(D_temp) > 0) {
      D_sample <- rbind(D_sample, D_temp[sample(1:nrow(D_temp), size = min(sample_size, nrow(D_temp)), replace = FALSE), ])
      qt1_mean <- mean(D_temp$qt1)
      qt2_mean <- mean(D_temp$qt2)
      Sigma <- stats::cov(D_temp[, c(1, 2)])
      Princip <- eigen(Sigma)
      flip_direction1 <- 0
      flip_direction2 <- 0
      if(Princip$vectors[1, 1] < 0){
        flip_direction1 <- 1
      }
      if(Princip$vectors[2, 2] < 0){
        flip_direction2 <- 1
      }
      M[,2 + 2 * i] <- qt1_mean + Princip$vectors[1, 1] * sqrt(Princip$values[1]) * cos(M$t) + Princip$vectors[1, 2] * sqrt(Princip$values[2]) * sin(M$t)
      M[,3 + 2 * i] <- qt2_mean + Princip$vectors[2, 1] * sqrt(Princip$values[1]) * cos(M$t) + Princip$vectors[2, 2] * sqrt(Princip$values[2]) * sin(M$t)
      Arrow_data[i + 1, 1] <- qt1_mean
      Arrow_data[i + 1, 2] <- qt2_mean
      Arrow_data[i + 1, 3] <- qt1_mean + (-1)^flip_direction1 * Princip$vectors[1,1] * sqrt(Princip$values[1])
      Arrow_data[i + 1, 4] <- qt2_mean + (-1)^flip_direction1 * Princip$vectors[2,1] * sqrt(Princip$values[1])
      Arrow_data[i + 4, 1] <- qt1_mean
      Arrow_data[i + 4, 2] <- qt2_mean
      Arrow_data[i + 4, 3] <- qt1_mean + (-1)^flip_direction2 * Princip$vectors[1,2] * sqrt(Princip$values[2])
      Arrow_data[i + 4, 4] <- qt2_mean + (-1)^flip_direction2 * Princip$vectors[2,2] * sqrt(Princip$values[2])
    }
  }
  ggplot2::ggplot(D_sample, ggplot2::aes(x = qt1 , y = qt2 ,color = g_factor))+
    ggplot2::geom_point()+ggplot2::theme_classic()+
    ggplot2::geom_smooth(method = 'lm', data = D, se = FALSE, formula = stats::as.formula('y ~ x')) +
    ggplot2::coord_fixed() +
    ggplot2::scale_color_manual(values = c('Non-carriers' = '#F8766D', 'Heterozygotes' = '#00BA38', 'Homozygotes' = '#619CFF')) +
    ggplot2::geom_segment(ggplot2::aes(x = Arrow_data[1, 1], y = Arrow_data[1, 2],
                              xend = Arrow_data[1, 3], yend = Arrow_data[1, 4] ),
                          color = 'red', size = 1, arrow = ggplot2::arrow()) +
    ggplot2::geom_segment(ggplot2::aes(x = Arrow_data[4, 1], y = Arrow_data[4, 2],
                              xend = Arrow_data[4, 3], yend = Arrow_data[4, 4] ),
                          color = 'red', size = 1, arrow = ggplot2::arrow()) +
    ggplot2::geom_segment(ggplot2::aes(x = Arrow_data[2, 1], y = Arrow_data[2, 2],
                              xend = Arrow_data[2, 3], yend = Arrow_data[2, 4] ),
                          color = 'green', size = 1, arrow = ggplot2::arrow()) +
    ggplot2::geom_segment(ggplot2::aes(x = Arrow_data[5, 1], y = Arrow_data[5, 2],
                              xend = Arrow_data[5, 3], yend = Arrow_data[5, 4] ),
                          color = 'green', size = 1, arrow = ggplot2::arrow()) +
    ggplot2::geom_segment(ggplot2::aes(x = Arrow_data[3, 1], y = Arrow_data[3, 2],
                              xend = Arrow_data[3, 3], yend = Arrow_data[3, 4] ),
                          color = 'blue', size = 1, arrow = ggplot2::arrow()) +
    ggplot2::geom_segment(ggplot2::aes(x = Arrow_data[6, 1], y = Arrow_data[6, 2],
                              xend = Arrow_data[6, 3], yend = Arrow_data[6, 4] ),
                          color = 'blue', size = 1, arrow = ggplot2::arrow()) +
    ggplot2::geom_polygon(data=M, ggplot2::aes(x=x0,y=y0), color='red',fill=NA,size=1.5) +
    ggplot2::geom_polygon(data=M, ggplot2::aes(x=x1,y=y1), color='green',fill=NA,size=1.5) +
    ggplot2::geom_polygon(data=M, ggplot2::aes(x=x2,y=y2), color='blue',fill=NA,size=1.5) +
    ggplot2::xlab(trait_name1)+ggplot2::ylab(trait_name2) + ggplot2::ggtitle(title)
}
