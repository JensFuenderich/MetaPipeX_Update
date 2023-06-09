#' Merging Replication Summaries
#'
#'
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#'
#' @description
#' \loadmathjax{}
#' \(
#' \\let\\underscore_
#' \)
#' Function to merge the replication statistics returned by MetaPipeX::create_replication_summaries() into a single data frame. This is the second function (and the fourth computational step) of the MetaPipeX pipeline. For more details on the pipeline, refer to the documentation of the MetaPipeX-package.
#'
#' @param data
#' The function expects the input to be a list of data frames or a path to a folder containing the replication summaries as .csv files. The input may either be produced by the MetaPipeX::create_replication_summaries() function, or any inputs that use the data template. A template of this data frame is available on \href{https://github.com/JensFuenderich/MetaPipeX/blob/main/Supplementary_Material/Table_Templates/2_Replication_Summaries/Replication_Summaries_template.csv}{{github}}, as is a \href{https://github.com/JensFuenderich/MetaPipeX/blob/main/Supplementary_Material/Table_Templates/2_Replication_Summaries/codebook_for_replication_summaries.csv}{{codebook}} for unambiguous identification of the abbreviations.
#'
#' @param output_folder
#' Define a path to which the merged replication summaries and the codebook are exported. If no path is specified, results are returned only in R.
#' @param suppress_list_output
#' A logical indicating whether results should be returned in R. If TRUE, no output is returned in R.
#'
#' @details
#' No transformations are performed on the data in this step of the MetaPipeX pipeline.
#'
#' @return
#' A list object containing the following components: \cr
#' ## merged_replication_summaries
#' A data frame with all replications from the input.
#'
#' ## codebook
#' A codebook that applies to the data frame (merged_replication_summaries). \cr
#' In order to export the data structure as .csv files in a folder, output_folder has to be specified.
#'
#' @examples
#'
#' # import the according table template
#' Replication_Summaries_template <- readr::read_csv(url(
#' paste("https://raw.githubusercontent.com/JensFuenderich/MetaPipeX/main/Supplementary_Material/",
#' "Table_Templates/2_Replication_Summaries/Replication_Summaries_template.csv",
#' sep = ""
#' )))
#'
#' # set seed for drawing data
#' set.seed(1973)
#'
#' # create vectors with names
#' MultiLab_names <- c("MultiLab_1", "MultiLab_1", "MultiLab_2",  "MultiLab_2")
#' ReplicationProject_names <- c("Effect_A", "Effect_B", "Effect_C", "Effect_D")
#' Replication_names <- c("Lab_A", "Lab_B", "Lab_C", "Lab_D", "Lab_E", "Lab_F", "Lab_G", "Lab_H")
#'
#' # random sampling for simulated data & building identifier variables
#' list_of_replication_summaries <- lapply(1:4, function(x){
#'   # sampling
#'   data_example <- as.data.frame(matrix(
#'   data = stats::rnorm(n = 200*(ncol(Replication_Summaries_template)-3), mean = 5, sd = 0.5),
#'   nrow = 200,
#'   ncol = ncol(Replication_Summaries_template)-3)
#'   )
#'   # rename columns according to template
#'   names(data_example) <- names(
#'   Replication_Summaries_template
#'   )[4:length(names(Replication_Summaries_template))]
#'   data_example$T_N <- round(data_example$T_N, 0)
#'   data_example$T_N <- round(data_example$C_N, 0)
#'   # building identifier variables
#'   MultiLab <- rep(MultiLab_names[x], times = nrow(data_example))
#'   ReplicationProject <- rep(ReplicationProject_names[x], times = nrow(data_example))
#'   Replication <- rep(if (x == 1 | x == 2) {
#'   Replication_names[1:4]
#'   } else if (x == 3 | x == 4) {
#'   Replication_names[5:8]
#'   }, each = nrow(data_example)/4)
#'   # combine data & identifiers
#'   cbind(MultiLab, ReplicationProject, Replication, data_example)
#' })
#' # rename list objects
#' names(list_of_replication_summaries) <- c("MultiLab_1_ReplicationProject_A_Replication_summaries",
#'                                           "MultiLab_1_ReplicationProject_B_Replication_summaries",
#'                                           "MultiLab_2_ReplicationProject_C_Replication_summaries",
#'                                           "MultiLab_2_ReplicationProject_D_Replication_summaries")
#'
#' ## applying the input to the MetaPipeX function
#'
#' # run merge_replication_summaries
#' example_MetaPipeX_output <- MetaPipeXUpdate::merge_replication_summaries(
#' data = list_of_replication_summaries
#' )
#'
#' \dontrun{
#' All examples with additional comments are available on github:
#' https://github.com/JensFuenderich/MetaPipeX/tree/main/Supplementary_Material/Code_Examples
#' }
#'
#'
#' @export
merge_replication_summaries <- function(data, output_folder = NULL, suppress_list_output = FALSE){

  ### Merge lab summaries

  if (is.list(data) == TRUE) {

    merged_replication_summaries <- do.call(rbind.data.frame, data)

  } else if(file.exists(data) == TRUE) {

    # collect file names from folder specified in data
    files <- list.files(path = file.path(data), pattern = "*.csv", full.names = T)

    # import the files and store as data frame
    tbl <- sapply(files, readr::read_csv, simplify=FALSE) %>% dplyr::bind_rows(.id = "id")
    merged_replication_summaries <- as.data.frame(tbl)

    # drop id column
    merged_replication_summaries$id <- NULL

    # drop redundant id column
    merged_replication_summaries$...1 <- NULL # not quite sure how to avoid this column being created in the first place

  } else {

    warning("Make sure to the 'data' input is either a list object or a valid path.")

  }

  ### Create codebook

  if (is.list(data) != TRUE && file.exists(data) != TRUE) {

    # function returns warning from the "Merge replication summaries" chunk, no further action necessary

  } else {

    # create empty df
    abbr_library <- data.frame(Abbreviation = logical(0),
                               Full_Name = logical(0))

    # pair abbreviations with verbal descriptions
    abbr_library <- as.data.frame(base::rbind(c("T_", "treatment group_"),
                                              c("C_", "control group_"),
                                              c("_N", "_number of participants"),
                                              c("_K", "_number of labs"),
                                              c("_MD", "_mean difference"),
                                              c("_Est_", "_model estimate for_"),
                                              c("_M", "_mean"),
                                              c("_SD", "_standard deviation"),
                                              c("SE_", "standard error of the_"),
                                              c("SMD", "standardized mean difference"),
                                              c("pooled_", "pooled_")
    ))

    # rename columns of df
    names(abbr_library) <- c("Abbreviation", "Full_Name")

    # extract names from merged df
    description_vector <- names(merged_replication_summaries)

    # sorry for this, did not want to loop
    # check if there's enough pipes in that orchestra
    #nrow(abbr_library) (the result of this should be equivalent to the max indexing in the following chunk)

    description_vector %<>%
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
      gsub(abbr_library$Abbreviation[11], abbr_library$Full_Name[11], .)

    description_vector <- gsub(pattern = "_", replacement = " ", description_vector)

    codebook_for_merged_replication_summeries <- data.frame(Variable_Name = names(merged_replication_summaries), Variable_Description = description_vector)
    codebook_for_merged_replication_summeries <- codebook_for_merged_replication_summeries[-c(1:3),]

    # do this one by hand, otherwise the abbr "MD" messes up the code
    codebook_for_merged_replication_summeries[codebook_for_merged_replication_summeries$Variable_Name == "MD",2] <- "mean difference"

    # add identifiers
    codebook_for_merged_replication_summeries <- rbind(data.frame(Variable_Name = c("MultiLab",
                                                                                    "ReplicationProject",
                                                                                    "Replication"),
                                                                  Variable_Description = c("The multi-lab in which the replication project was publicised (e.g., ML2)",
                                                                                           "The name of the replication project (or replicated target-effect)",
                                                                                           "The replication (e.g., lab name) that a data point is associated with")),
                                                       codebook_for_merged_replication_summeries)

  }

  ## Outputs

  if (is.null(output_folder) == TRUE) {

    base::print("You chose not to export the data as .csv files.")

  } else {

    # export .csv files
    readr::write_csv(merged_replication_summaries,
                     paste(output_folder, "Merged_Replication_Summaries.csv", sep = ""))
    readr::write_csv(codebook_for_merged_replication_summeries,
                     paste(output_folder, "codebook_for_merged_replication_summeries.csv", sep = ""))

  }

  if (suppress_list_output == TRUE) {

    base::print("You chose not to return results in R. If you specified an output folder, check that folder for the code book and merged replication summaries.")

  } else if (suppress_list_output == FALSE) {

    # create list output
    output <- list(merged_replication_summaries, codebook_for_merged_replication_summeries)

    # rename list elements
    names(output) <- c("Merged_Replication_Summaries", "codebook_for_merged_replication_summeries")

    # return the output (function aborts here)
    return(output)

  }

}
