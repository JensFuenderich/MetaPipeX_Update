#' Meta Analyses per MASC (meta-analytical study collection)
#'
#' @import mathjaxr
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#'
#' @description
#' \loadmathjax{}
#' \(
#' \\let\\underscore_
#' \)
#' Function to run meta-analyses on the mean difference (MD) and the standardized mean difference (SMD). The meta-analyses are run with the metafor::rma.mv function (Viechtbauer, 2010). For more details on the meta-analyses, refer to the Details and Return section. This function is the third (and fifth computational step) of the MetaPipeX pipeline. For more details on the pipeline, refer to the documentation of the MetaPipeX-package.
#'
#' @param data The function expects the input to be a data frame. The input may either be the data frame produced by the MetaPipeXUpdate::merge_site_summaries() function, or one with the same columns names. A template of this data frame is available on \href{https://github.com/JensFuenderich/MetaPipeX/blob/main/Supplementary_Material/Table_Templates/3_Merged_Site_Summaries/Merged_Site_Summaries_template.csv}{{github}}, as is a \href{https://github.com/JensFuenderich/MetaPipeX/blob/main/Supplementary_Material/Table_Templates/3_Merged_Site_Summaries/codebook_for_merged_site_summeries.csv}{{codebook}} for unambiguous identification of the abbreviations. Further, it is possible to use a \href{INSERT LINK}{{reduced version of the codebook}}, as meta-analyses are applied to MD and SMD only.
#' @param output_folder Specify the output folder for the meta-analyses and the codebook. If no folder is specified, the function will return its output only to the R environment (unless this is suppressed under suppress_list_output).
#' @param suppress_list_output A logical indicating whether results should be returned in R. If TRUE, no output is returned in R.
#' @param method A character string to specify the type of model to be fitted. Default is “REML”. For more details, refer to the \href{https://www.metafor-project.org/doku.php/help}{{metafor}}  documentation.
#' @param sparse A logical indicating whether sparse matrices should be used.
#'
#' @details
#'
#' The meta-analyses within the function are written with metafor::rma.mv (Viechtbauer, 2010). The multivariate version of the rma function is deployed to allow for the use of sparse matrices (“sparse = TRUE”) for optimal performance in meta-analyses with thousands of data collection sites. They are fitted as a random-effects model with “random = ~ 1 | Data_Collection_Site” and a restricted maximum likelihood estimation (“REML”).
#' The function runs seven meta-analyses per MASC:
#' \itemize{
#'  \item{control mean (yi = C_M, V = SE_C_M)}
#'  \item{treatment mean (yi = T_M, V = SE_T_M)}
#'  \item{control group standard deviation(yi = ln_SD, V = SE_ln_SD)}
#'  \item{treatment group standard deviation (yi = ln_SD, V = SE_ln_SD)}
#'  \item{mean difference (yi = MD, V = SE_MD)}
#'  \item{pooled standard deviation (yi = ln_SD, V = SE_ln_SD)}
#'  \item{standardized mean difference (yi = SMD, V = SE_SMD)}
#'  }
#'
#'  Meta-analyses of standard deviations are performed on log-transformed standard deviations, as recommended by Nakagawa et al. (2015). All standard deviations (C_SD, T_SD, pooled_SD) are transformed using:
#'  \itemize{
#'  \item{R-Code} \cr
#'    \code{ln_SD = log(SD) + 1 / (2 * (n - 1)) } \cr
#'  \item{Model} \cr
#'    {
#'      \mjdeqn{ ln \hat{\tau} =  ln SD + \frac{1}{2(n-1)} }{}
#'    }
#'  }
#'  All standard errors of standard deviations are created using:
#' \itemize{
#'  \item{R-Code} \cr
#'    \code{SE_ln_SD = sqrt(1 / (2 * (n - 1))) } \cr
#'  \item{Model} \cr
#'    {
#'      \mjdeqn{  s_{ln \hat{\tau}}^2 = \sqrt{ \frac{1}{2(n-1)} } }{}
#'    }
#'  }
#'
#'  In order to achieve interpretable meta-analytical results (the model estimate, tau, tau^2 and the coefficient of variation), they are transformed back into the original units. The following transformation is applied to the model estimate:
#'  \itemize{
#'  \item{R-Code} \cr
#'    \code{Model_Est = exp(mu_ln + 0.5 * tau2_ln)} \cr
#'  \item{Model} \cr
#'    {
#'      \mjdeqn{  \mu = e^{\mu_{ln}+0.5\tau_{ln}^2} }{}
#'    }
#'  }
#'  The following transformation is applied to tau^2 (which is subsequently used to calculate tau and the coefficient of variation):
#'  \itemize{
#'  \item{R-Code} \cr
#'    \code{Tau2 =  exp(2 * mu_ln + tau2_ln) * (exp(tau2_ln) - 1) } \cr
#'  \item{Model} \cr
#'    {
#'      \mjdeqn{  \tau^2 = e^{2 \mu_{ln} + \tau_{ln}^2} ( e^{\tau_{ln}^2} - 1 ) }{}
#'    }
#'  }
#'
#'
#'
#' @return
#' The output is a list with two objects: A data frame with the meta-analytical results and a codebook for unambiguous identification of its columns. \cr
#' \cr ## meta analyses \cr
#' The data frame contains information to identify each analysis (MultiLab, MASC) and statistical output from the two meta-analyses per MASC. The statistical output for each meta-analysis includes:
#' \itemize{
#' \item{A model estimate for the y of interest (Est__).}
#' \item{The number of sites included in the analysis (Result__K).}
#' \item{The estimated \mjeqn{\tau^2}{} (sigma2 from the rma.mv object) value (Tau2__).}
#' \item{The estimated \mjeqn{\tau}{} (the square root of the sigma2 from the rma.mv object) value (Tau2__).}
#' \item{The estimated \mjeqn{I^2}{} value. \mjeqn{I^2}{} is not part of the rma.mv output object and has to be calculated from \mjeqn{\tau}{}.
#' \mjdeqn{ I^2 = 100 \frac{ \hat{\tau}^2  }{ \hat{\tau}^2 + \tilde{v}} }{}
#' with
#' \mjdeqn{ \tilde{v} = \frac{(k-1)\sum w_{i}}{\left(\sum  w_{i}\right)-\sum w_{i}^2} }{}
#' Transformation according to: https://wviechtb.github.io/metafor/reference/print.rma.html
#' }
#' \item{The estimated \mjeqn{H^2}{} value. \mjeqn{H^2}{} is not part of the rma.mv output object and has to be calculated from \mjeqn{\tau}{}.
#' \mjdeqn{ H^2 = 100 \frac{ \hat{\tau}^2 + \tilde{v} }{\tilde{v}} }{}
#' with
#' \mjdeqn{ \tilde{v} = \frac{(k-1)\sum w_{i}}{\left(\sum  w_{i}\right)-\sum w_{i}^2} }{}
#' }
#' \item{The Q statistic (QE__).}
#' \item{The p-value from the test on the Q statistic (QEp__).}
#' }
#' ## codebook \cr
#' A codebook that applies to the data frame (meta_analyses).
#'
#'
#' @references
#'
#' Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. Journal of Statistical Software, 36(3), 1-48. doi: 10.18637/jss.v036.i03
#' Nakagawa, S., Poulin, R., Mengersen, K., Reinhold, K., Engqvist, L., Lagisz, M., & Senior, A. M. (2015). Meta‐analysis of variation: ecological and evolutionary applications and beyond. Methods in Ecology and Evolution, 6(2), 143-152. doi: 10.1111/2041-210X.12309
#'
#' @examples
#'
#' # create IPD for 10 MASCs (all from the same MultiLab)
#' sim_out <- mapply(MetaPipeXUpdate::simulate_IPD,
#'                   MASC_index = 1:5,
#'                   seed = 50 + (0:(5-1)),
#'                   SIMPLIFY = FALSE)
#' # rename list elements (the individual MASCs)
#' names(sim_out) <- paste("MASC", 1:5, sep = "")
#'
#' # create site summaries
#' Site_Summaries <- MetaPipeXUpdate::summarize_sites(data = sim_out)
#'
#' # merge site summaries
#' Merged_Site_Summaries <- MetaPipeXUpdate::merge_site_summaries(data = Site_Summaries$Site_Summaries)
#'
#' # run the MetaPipeX function to meta-analyze all MASCs
#' MetaPipeXUpdate::meta_analyze_MASCs(data = Merged_Site_Summaries$Merged_Site_Summaries)
#'
#' \dontrun{
#' All examples with additional comments are available on github:
#' https://github.com/JensFuenderich/MetaPipeX/tree/main/Supplementary_Material/Code_Examples
#' }
#'
#' @export
meta_analyze_MASCs <- function(data, output_folder = NULL, suppress_list_output = FALSE, method = "REML", sparse = FALSE){

  ## input is a large df with all multi-labs & MASCs


  ### Run meta-analyses

  ## create a function that runs all meta-analyses for one MASC (e.g., target-effect)

  single_MASC_analyses <- function(subset_MASC){

    method <- method

    # create a vector with the column names for the analysis
    col_names <- c(
      # N per MASC & K (Number of Sites):
      "Result__K",
      "Result__N",
      # Meta-analytic Estimates:
      "Est__C_M",
      "Est__T_M",
      "Est__C_SD",
      "Est__T_SD",
      "Est__pooled_SD",
      "Est__MD",
      "Est__SMD",
      # pval of Meta-analytic Estimates:
      "pval_Est__C_M",
      "pval_Est__T_M",
      "pval_Est__C_SD",
      "pval_Est__T_SD",
      "pval_Est__pooled_SD",
      "pval_Est__MD",
      "pval_Est__SMD",
      # K of Meta-analytic Estimates:
      "Est__C_M_K",
      "Est__T_M_K",
      "Est__C_SD_K",
      "Est__T_SD_K",
      "Est__pooled_SD_K",
      "Est__MD_K",
      "Est__SMD_K",
      # Tau2 Analyses and SE of Tau2:
      "Tau2__C_M",
      "Tau2__T_M",
      "Tau2__C_SD",
      "Tau2__T_SD",
      "Tau2__pooled_SD",
      "Tau2__MD",
      "Tau2__SMD",
      # Tau Analyses:
      "Tau__C_M",
      "Tau__T_M",
      "Tau__C_SD",
      "Tau__T_SD",
      "Tau__pooled_SD",
      "Tau__MD",
      "Tau__SMD",
      # Coefficient of Variation Analyses:
      "CoeffVar__T_M",
      "CoeffVar__C_M",
      "CoeffVar__T_SD",
      "CoeffVar__C_SD",
      "CoeffVar__pooled_SD",
      "CoeffVar__MD",
      "CoeffVar__SMD",
      # I2 Analyses:
      "I2__T_M",
      "I2__C_M",
      "I2__T_SD",
      "I2__C_SD",
      "I2__MD",
      "I2__pooled_SD",
      "I2__SMD",
      # H2 Analyses:
      "H2__T_M",
      "H2__C_M",
      "H2__T_SD",
      "H2__C_SD",
      "H2__MD",
      "H2__pooled_SD",
      "H2__SMD",
      # QE Values:
      "QE__T_M",
      "QE__C_M",
      "QE__T_SD",
      "QE__C_SD",
      "QE__MD",
      "QE__pooled_SD",
      "QE__SMD",
      # QEp Values:
      "QEp__T_M",
      "QEp__C_M",
      "QEp__T_SD",
      "QEp__C_SD",
      "QEp__MD",
      "QEp__pooled_SD",
      "QEp__SMD"
    )


    # create a df for the results of the analysis from the subset
    MASC.df <- data.frame(t(rep(0,length(col_names))))
    # rename columns
    names(MASC.df) <- col_names

    ## replace infinite values (Inf) in input df with NA
    subset_MASC <- do.call(data.frame, lapply(subset_MASC, function(value){replace(value, is.infinite(value),NA)}))

    ## insert information on sample sizes and number of data collection sites

    # N
    MASC.df["Result__N"] <- sum(subset_MASC$T_N + subset_MASC$C_N)
    # K
    MASC.df["Result__K"] <- length(subset_MASC$Data_Collection_Site)

    ## Transformations for rma.mv output, which does not include I2 and H2
    # transformations according to https://cran.r-project.org/web/packages/metafor/metafor.pdf
    # I2 as described in https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate

    I2_fct <- function(rma_mv_obj){
      k <- rma_mv_obj$k
      wi <- 1/rma_mv_obj$vi
      vt <- (k-1) * sum(wi) / (sum(wi)^2 - sum(wi^2))
      100 * rma_mv_obj$sigma2 / (rma_mv_obj$sigma2 + vt)
    }

    H2_fct <- function(rma_mv_obj){
      k <- rma_mv_obj$k
      wi <- 1/rma_mv_obj$vi
      vt <- (k-1) * sum(wi) / (sum(wi)^2 - sum(wi^2))
      (rma_mv_obj$sigma2 + vt)/vt
    }

    ## log transformation for meta-analyses of standard deviations
    # transformations according to Nakagawa et al. (2015) "Meta‐analysis of variation")

    # define ln of standard deviation
    ln_SD_fun <- function(SD, n) {
      base::log(SD) + 1 / (2 * (n - 1))
    }
    # define standard error for ln of standard deviation
    SE_lnSD_fun <- function(n) {
      sqrt(1 / (2 * (n - 1)))
    }

    ## transformation back from log to original unit for meta-analyses of standard deviations
    # transformations according to https://mathworld.wolfram.com/LogNormalDistribution.html

    # define transformation to original unit for the model estimate
    Model_Est_fun <- function(mu_ln, tau2_ln){
      exp(mu_ln + 0.5 * tau2_ln)
    }
    # define transformation to original unit for tau2
    Tau2_fun <- function(mu_ln, tau2_ln){
      exp(2 * mu_ln + tau2_ln)* (exp(tau2_ln) - 1)
    }

    ## run meta-analyses and fill "MASC.df" with the output

    # 1 Heterogeneity of control mean
    MASC.df <- MetaPipeXUpdate:::standard_MA_subset(MASC_data = subset_MASC,
                                                    Data_Collection_Site = "Data_Collection_Site",
                                                    yi = "C_M",
                                                    SE = "SE_C_M",
                                                    method = method,
                                                    sparse = sparse,
                                                    MASC.df = MASC.df)

    # 2 Heterogeneity of treatment mean
    MASC.df <- MetaPipeXUpdate:::standard_MA_subset(MASC_data = subset_MASC,
                                                    Data_Collection_Site = "Data_Collection_Site",
                                                    yi = "T_M",
                                                    SE = "SE_T_M",
                                                    method = method,
                                                    sparse = sparse,
                                                    MASC.df = MASC.df)

    # 3 Heterogeneity of control group SD
    MASC.df <- MetaPipeXUpdate:::ln_MA_subset(MASC_data = subset_MASC,
                                              Data_Collection_Site = "Data_Collection_Site",
                                              yi = "C_SD",
                                              SE = "SE_C_SD",
                                              N = "C_N",
                                              method = method,
                                              sparse = sparse,
                                              MASC.df = MASC.df)


    # 4 Heterogeneity of treatment group SD
    MASC.df <- MetaPipeXUpdate:::ln_MA_subset(MASC_data = subset_MASC,
                                              Data_Collection_Site = "Data_Collection_Site",
                                              yi = "T_SD",
                                              SE = "SE_T_SD",
                                              N = "T_N",
                                              method = method,
                                              sparse = sparse,
                                              MASC.df = MASC.df)


    # 5 Heterogeneity of mean difference
    MASC.df <- MetaPipeXUpdate:::standard_MA_subset(MASC_data = subset_MASC,
                                                    Data_Collection_Site = "Data_Collection_Site",
                                                    yi = "MD",
                                                    SE = "SE_MD",
                                                    method = method,
                                                    sparse = sparse,
                                                    MASC.df = MASC.df)


    # 6 Heterogeneity of pooled SD
    MASC.df <- MetaPipeXUpdate:::ln_MA_subset(MASC_data = subset_MASC,
                                              Data_Collection_Site = "Data_Collection_Site",
                                              yi = "pooled_SD",
                                              SE = "SE_pooled_SD",
                                              N = "C_N",
                                              N_second_group = "T_N",
                                              method = method,
                                              sparse = sparse,
                                              MASC.df = MASC.df)


    # 7 Heterogeneity of effect size g (Borenstein)
    MASC.df <- MetaPipeXUpdate:::standard_MA_subset(MASC_data = subset_MASC,
                                                    Data_Collection_Site = "Data_Collection_Site",
                                                    yi = "SMD",
                                                    SE = "SE_SMD",
                                                    method = method,
                                                    sparse = sparse,
                                                    MASC.df = MASC.df)

    # add descriptive columns
    descriptive_columns <- data.frame(MultiLab = unique(subset_MASC$MultiLab),
                                      MASC = unique(subset_MASC$MASC))
    MASC.df <- cbind(descriptive_columns, MASC.df)

    return(MASC.df)

  }

  ## prepare data for analyses

  # create list: 1 list object = 1 multilab
  data_list_MultiLab_split <- split(data, data$MultiLab)

  # create nested list with the MASCs as level 2 list objects & rename list output
  nested_data_list_MASC_split <- lapply(data_list_MultiLab_split, function(x){split(x, x$MASC)})
  names(nested_data_list_MASC_split) <- names(data_list_MultiLab_split)

  # apply the meta-analyses & rename list output
  nested_list_output <- lapply(nested_data_list_MASC_split,function(x){lapply(x, single_MASC_analyses)})
  names(nested_list_output) <- names(data_list_MultiLab_split)

  # create output df from list object
  meta_analyses <- dplyr::bind_rows(lapply(nested_list_output, function(x){dplyr::bind_rows(x)}))

  ### Create codebook

  # create empty df
  abbr_library <- data.frame(Abbreviation = logical(0),
                             Full_Name = logical(0))

  # pair abbreviations with verbal descriptions
  abbr_library <- as.data.frame(base::rbind(c("_N", "_number of participants"),
                                            c("_K", "_number of data collection sites"),
                                            c("_C_", "__control group_"),
                                            c("_T_", "__treatment group_"),
                                            c("__pooled_", "__pooled_"),
                                            c("_M", "_mean"),
                                            c("_SD", "_standard deviation"),
                                            c("_MD", "_mean difference"),
                                            c("_SMD", "_standardized mean difference"),
                                            c("__SE_", "__standard error of the_"),
                                            c("_Est_", "_model estimate for_"),
                                            c("Tau2__", "Tau2 for__"),
                                            c("Tau__", "Tau for__"),
                                            c("CoeffVar__", "Coefficient of Variation (tau/mu) for__"),
                                            c("I2__", "I2 for__"),
                                            c("H2__", "H2 for__"),
                                            c("QE__", "QE for__"),
                                            c("QEp__", "QEp for__")
  ))

  # rename columns of df
  names(abbr_library) <- c("Abbreviation", "Full_Name")

  # extract names from merged df
  description_vector <- names(meta_analyses)

  # sorry for this, did not want to loop
  # check if there's enough pipes in that orchestra
  #nrow(abbr_library) (the result of this should be equivalent to the max indexing in the following chunk)


  description_vector %<>% # pipe from magrittr
    gsub(abbr_library$Abbreviation[1], abbr_library$Full_Name[1], .) %>%
    gsub(abbr_library$Abbreviation[2], abbr_library$Full_Name[2], .) %>%
    gsub(abbr_library$Abbreviation[3], abbr_library$Full_Name[3], .) %>%
    gsub(abbr_library$Abbreviation[4], abbr_library$Full_Name[4], .) %>%
    gsub(abbr_library$Abbreviation[5], abbr_library$Full_Name[5], .) %>%
    gsub(abbr_library$Abbreviation[6], abbr_library$Full_Name[6], .) %>%
    gsub(abbr_library$Abbreviation[7], abbr_library$Full_Name[7], .) %>%
    gsub(abbr_library$Abbreviation[8], abbr_library$Full_Name[8], .) %>%
    gsub(abbr_library$Abbreviation[9], abbr_library$Full_Name[9], .) %>%
    gsub(abbr_library$Abbreviation[10], abbr_library$Full_Name[10], .) %>%
    gsub(abbr_library$Abbreviation[11], abbr_library$Full_Name[11], .) %>%
    gsub(abbr_library$Abbreviation[12], abbr_library$Full_Name[12], .) %>%
    gsub(abbr_library$Abbreviation[13], abbr_library$Full_Name[13], .) %>%
    gsub(abbr_library$Abbreviation[14], abbr_library$Full_Name[14], .) %>%
    gsub(abbr_library$Abbreviation[15], abbr_library$Full_Name[15], .) %>%
    gsub(abbr_library$Abbreviation[16], abbr_library$Full_Name[16], .) %>%
    gsub(abbr_library$Abbreviation[17], abbr_library$Full_Name[17], .) %>%
    gsub(abbr_library$Abbreviation[18], abbr_library$Full_Name[18], .)

  description_vector <- sub(pattern = "__Result__", replacement = "_", description_vector)
  description_vector <- sub(pattern = "___", replacement = "_", description_vector)
  description_vector <- sub(pattern = "__", replacement = "_", description_vector)
  description_vector <- sub(pattern = "_", replacement = " ", description_vector)

  codebook_for_meta_analyses <- data.frame(Variable_Name = names(meta_analyses), Variable_Description = description_vector)

  ## Outputs

  if (is.null(output_folder)) {

    base::print("You chose not to export the data as .csv files.")

  } else {

    # export .csv files
    readr::write_csv(meta_analyses,
                     paste(output_folder, "meta_analyses.csv", sep = ""))
    readr::write_csv(codebook_for_meta_analyses,
                     paste(output_folder, "codebook_for_meta_analyses.csv", sep = ""))

  }

  if (suppress_list_output == TRUE) {

    base::print("You chose not to return results in R. If you specified an output folder, check that folder for the code book and meta-analyses output.")

  } else if (suppress_list_output == FALSE) {

    # create list output
    output <- list(meta_analyses, codebook_for_meta_analyses)

    # rename list elements
    names(output) <- c("Meta_Analyses", "codebook_for_meta_analyses")

    # return the output (function aborts here)
    return(output)

  }


}
