## log transformation for meta-analyses of standard deviations
# transformations according to Nakagawa et al. (2015) "Meta‐analysis of variation")

# define standard error for ln of standard deviation
SE_lnSD_fun <- function(n) {
  sqrt(1 / (2 * (n - 1)))
}
