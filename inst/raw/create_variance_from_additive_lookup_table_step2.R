#############################################################################
# Collect files from step 1 to create variance lookup object
#############################################################################
library(data.table)

# Note that these are exactly the same values that were defined in step1
f_alleles <- c(seq(0.001, 0.01, by = 0.001), seq(0.02,0.1, by = 0.01), seq(0.15,0.5, by = 0.05))

# Read all the results from step1
.variance.simulation <- rbindlist(lapply(1:length(f_alleles), function(x){
  out_dat <- fread(paste0("./job_", x, ".txt"))
  return(out_dat)
}))

# We save this in a place for package internal data
usethis::use_data(.variance.simulation, internal = TRUE, overwrite = TRUE)

library(ggplot2)
library(ggsci)

# One way to visualise this
g_p <- ggplot(.variance.simulation, aes(x = beta_hat, y = alpha_hat, color = factor(f_par), group = factor(f_par))) +
  geom_line() + theme_bw() + scale_color_viridis_d() + xlab("Estimated beta after inv norm") + ylab("Estimated alpha after inv norm") +
  theme(legend.position = "top") + labs(color="allele frequency") +
  ggtitle("Beta hat and Alpha hat")
