## log transformation for meta-analyses of standard deviations
# transformations according to Nakagawa et al. (2015) "Meta‐analysis of variation")

# define ln of standard deviation
ln_SD_fun <- function(SD, n) {
  base::log(SD) + 1 / (2 * (n - 1))
}
