## function for the pipeline: runs a single meta-analysis for standard deviations using metafor::rma.mv
ln_MA_subset <- function(MASC_data, Data_Collection_Site, yi, SE, N, N_second_group = NULL, method, sparse, MASC.df){

  if ( nrow(stats::na.omit(MASC_data[, c(Data_Collection_Site, yi, SE, N, N_second_group)])) <= 1 ) {} else {

    if (is.null(N_second_group)) {
      # use the subset of columns relevant to this analysis
      MASC_data <- stats::na.omit(MASC_data[, c(Data_Collection_Site, yi, SE, N)])
      names(MASC_data) <- c("Data_Collection_Site", "yi", "SE", "N")
      # create df with ln data
      ln_data_full <- data.frame(
        Data_Collection_Site = MASC_data$Data_Collection_Site,
        ln_SD = MetaPipeXUpdate:::ln_SD_fun(SD = MASC_data$yi,
                                            n = MASC_data$N),
        SE_ln_SD = MetaPipeXUpdate:::SE_lnSD_fun(n = MASC_data$N)
      )
    } else {
      # use the subset of columns relevant to this analysis
      MASC_data <- stats::na.omit(MASC_data[, c(Data_Collection_Site, yi, SE, N, N_second_group)])
      names(MASC_data) <- c("Data_Collection_Site", "yi", "SE", "N", "N_second_group")
      # create df with ln data
      ln_data_full <- data.frame(
        Data_Collection_Site = MASC_data$Data_Collection_Site,
        ln_SD = MetaPipeXUpdate:::ln_SD_fun(SD = MASC_data$yi,
                                            n = MASC_data$N + MASC_data$N_second_group),
        SE_ln_SD = MetaPipeXUpdate:::SE_lnSD_fun(n = MASC_data$N + MASC_data$N_second_group)
      )
    }

    # remove values that irritate the meta-analysis
    ln_data <- ln_data_full %>% dplyr::filter(is.finite(ln_SD))
    # print warning, if necessary
    if (nrow(ln_data_full) > nrow(ln_data)) {
      print("data collection sites were removed, due to non-positive ln(C_SD)")
    }
    # run the meta-analysis
    rma_output <- metafor::rma.mv(yi = ln_SD,
                                  V = SE_ln_SD^2,
                                  random = ~ 1 | Data_Collection_Site,
                                  method = method,
                                  sparse = TRUE,
                                  data = ln_data)
    # insert the meta analysical results at the appropriate columns in the df
    # transformations for Est_, Tau_ and /Tau2_ according to:
    # https://stats.stackexchange.com/questions/241187/calculating-standard-deviation-after-log-transformation
    # insert the meta analysical results at the appropriate columns in the df
    MASC.df[paste0("Est__", yi)] <- MetaPipeXUpdate:::Model_Est_fun(as.numeric(rma_output$b), rma_output$sigma2)
    MASC.df[paste0("Est__", yi, "_K")] <- rma_output$k
    MASC.df[paste0("pval_Est__", yi)] <- rma_output$pval
    MASC.df[paste0("Tau2__", yi)] <- MetaPipeXUpdate:::Tau2_fun(as.numeric(rma_output$b), rma_output$sigma2)
    MASC.df[paste0("Tau__", yi)] <- sqrt(MASC.df[paste0("Tau2__", yi)])
    MASC.df[paste0("CoeffVar__", yi)] <- sqrt(MetaPipeXUpdate:::Tau2_fun(as.numeric(rma_output$b), rma_output$sigma2))/abs(MetaPipeXUpdate:::Model_Est_fun(as.numeric(rma_output$b), rma_output$sigma2))
    MASC.df[paste0("I2__", yi)] <- MetaPipeXUpdate:::I2_fct(rma_mv_obj = rma_output)
    MASC.df[paste0("H2__", yi)] <- MetaPipeXUpdate:::H2_fct(rma_mv_obj = rma_output)
    MASC.df[paste0("QE__", yi)] <- rma_output$QE
    MASC.df[paste0("QEp__", yi)] <- rma_output$QEp

    rm(MASC_data, ln_data, rma_output)

    MASC.df
  }

}
