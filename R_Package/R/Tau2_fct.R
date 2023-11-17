## transformation back from log to original unit for meta-analyses of standard deviations
# transformations according to https://mathworld.wolfram.com/LogNormalDistribution.html

# define transformation to original unit for tau2
Tau2_fun <- function(mu_ln, tau2_ln){
  exp(2 * mu_ln + tau2_ln)* (exp(tau2_ln) - 1)
}
