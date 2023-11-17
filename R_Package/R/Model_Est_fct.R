## transformation back from log to original unit for meta-analyses of standard deviations
# transformations according to https://mathworld.wolfram.com/LogNormalDistribution.html

# define transformation to original unit for the model estimate
Model_Est_fun <- function(mu_ln, tau2_ln){
  exp(mu_ln + 0.5 * tau2_ln)
}
