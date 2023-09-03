#############################################################################
# Variance effects can be created from strong additive effects
# We can simulate this and interpolate this data to quickly assess the
# potential impact for a given MAF/Beta combination
# To run this, call the script with:
# JOB_IND=1 Rscript ./inst/raw/create_variance_from_additive_lookup_table_step1.R
# For JOB_IND=1 up to JOB_IND=27
# This will create the 27 job files that are parsed in step2
#############################################################################
library(data.table)
library(ggplot2)

# Simple inverse normal transform
std <- function(x) {
  I <- !is.na(x)
  r <- rank(x[I])
  z <- qnorm(r/(sum(I) + 1))
  y <- x
  y[I] <- z
  return(y)
}


# This is currently implemented in the package, but this is the version that
# was used for this simulation
alpha.calc <- function(qt,g){
  n <- length(qt)
  f <- mean(g) / 2.0

  qt0 <- qt[g == 0]
  qt1 <- qt[g == 1]
  qt2 <- qt[g == 2]
  m0 <- mean(qt0)
  m1 <- mean(qt1)
  m2 <- mean(qt2)
  rss0 <- sum((qt0 - m0)^2)
  rss1 <- sum((qt1 - m1)^2)
  rss2 <- sum((qt2 - m2)^2)

  # Estimators for null-model and log-likelihood function
  sigma2_null <- (rss0 + rss1 + rss2) / n
  l0 <- -n * log(sigma2_null) / 2.0

  # Estimators for alt-model and log-likelihood function
  a <- 2*f*rss0
  b <- (2*f-1)*rss1
  c <- (2*f-2)*rss2

  alpha <- (-b+sqrt(b^2-4*a*c))/(2*a)
  sigma2_alt <- (rss2/(alpha^2)+rss1/alpha + rss0)/n

  l1 <- -n*log(sigma2_alt)/2-log(alpha)*n*f

  # Calculation of significance
  X2 <- 2*(l1-l0)
  pval <- pchisq(X2,1,lower.tail=F)
  return(list(alpha = alpha, sigma2_alt = sigma2_alt, pval = pval))
}

# Grid of frequencies for the simulation
f_alleles <- c(seq(0.001, 0.01, by = 0.001), seq(0.02,0.1, by = 0.01), seq(0.15,0.5, by = 0.05))

# We run each of these in a specific job, since this is a bit heavy
job_ind <- as.numeric(Sys.getenv(c("JOB_IND")))
f_allele <- f_alleles[job_ind]

# Simulate the variance effect that is expected from an additive effect
betas <- seq(0,5, length.out = 500)
samp_size <- 100000 # High enough to get enough for rare...
n_sim <- 500

beta_hat_vec <- rep(0, length(betas))
alpha_hat_vec <- rep(0, length(betas))
alpha_hat_std_vec <- rep(0, length(betas))
for(i in 1:length(betas)){
  print(i)
  # Init
  beta_hat_sim <- rep(0, n_sim)
  alpha_hat_sim <- rep(0, n_sim)
  alpha_hat_std_sim <- rep(0, n_sim)
  for(j in 1:n_sim){
    genotype1 <- sample(c(0,1),
                        size = samp_size,
                        replace = TRUE,
                        prob = c(1-f_allele,f_allele))
    genotype2 <- sample(c(0,1),
                        size = samp_size,
                        replace = TRUE,
                        prob = c(1-f_allele,f_allele))
    genotype <- genotype1 + genotype2
    my_samp <- rnorm(samp_size) +
      betas[i]*genotype
    my_samp_int <- dcutils::std(my_samp)
    fm <- lm(my_samp_int~genotype)
    my_samp_corr_int <- std(residuals(fm))
    beta_hat_sim[j] <- summary(fm)$coefficients[2]
    alpha_hat_sim[j] <- alpha.calc(my_samp_int, genotype)$alpha
    alpha_hat_std_sim[j] <- alpha.calc(my_samp_corr_int, genotype)$alpha
  }
  beta_hat_vec[i] <- mean(beta_hat_sim)
  alpha_hat_vec[i] <- mean(alpha_hat_sim)
  alpha_hat_std_vec[i] <- mean(alpha_hat_std_sim)
}

out_df <- data.frame(f_par = f_allele,
                     beta_hat = beta_hat_vec,
                     alpha_hat = alpha_hat_vec,
                     alpha_hat_std = alpha_hat_std_vec,
                     beta_orig = betas)

# This only needs to be done once, and CRAN does not allow any writing
# into userspace, so we comment this line out for now.
#data.table::fwrite(out_df, file = paste0("./job_", job_ind, ".txt"))
