% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_replication_summaries.R
\name{create_replication_summaries}
\alias{create_replication_summaries}
\title{Creating Replication Summaries}
\usage{
create_replication_summaries(
  data,
  MultiLab = NULL,
  ReplicationProject = NULL,
  Replication = NULL,
  DV = NULL,
  Group = NULL,
  output_folder = NULL,
  suppress_list_output = FALSE
)
}
\arguments{
\item{data}{A data frame or list of data frames that contain the individual participant data. The function expects the relevant columns to be named consistently across all list objects. Relevant to this function are columns that represent information on the MultiLab (e.g., Many Labs 2), the ReplicationProject (e.g., Ross1), the Replication (the lab a data point is assigned to), the group (either the treatment or control condition) and the single data point of the dependent variable (DV) per person. A template of this data frame is available on \href{https://github.com/JensFuenderich/MetaPipeX/blob/main/Supplementary_Material/Table_Templates/1_Individual_Participant_Data/IPD_template.csv}{{github}}, as is a \href{https://github.com/JensFuenderich/MetaPipeX/blob/main/Supplementary_Material/Table_Templates/1_Individual_Participant_Data/codebook_for_individual_participant_data.csv}{{codebook}} for unambiguous identification of the abbreviations.}

\item{MultiLab}{Character vector with the name of the columns in the list elements of "data" that contain the multi-lab name(s). If \emph{is.null(MultiLab) == TRUE}, "MultiLab" is chosen as the default.}

\item{ReplicationProject}{Character vector with the name of the columns in the list elements of "data" that contain the replication project name(s). If \emph{is.null(ReplicationProject) == TRUE}, "ReplicationProject" is chosen as the default. Each replication project comprises a single target-effect with direct replications across multiple replications (/labs).}

\item{Replication}{Character vector with the name of the columns in the list elements of "data" that contain the replication names. If \emph{is.null(Replication) == TRUE}, "Replication" is chosen as the default. The meta-analyses in MetaPipeX::meta_analyses() and MetaPipeX::full_pipeline() are run as random effects models in metafor::rma.mv() with “random = ~ 1 | Replication”. Thus, the pipeline assumes a distribution of true statistics (e.g., treatment means, mean differences, standardized mean differences).}

\item{DV}{Character vector with the name of the columns in the list elements of "data" that contain the (aggregated) dependent variable. If \emph{is.null(DV) == TRUE}, "DV" is chosen as the default.}

\item{Group}{Character vector with the name of the columns in the list elements of "data" that contain the (treatment/control) group identification. If \emph{is.null(Group) == TRUE}, "Group" is chosen as the default. These should only contain values of 0 (control group), 1 (treatment group) and NA (unidentified).}

\item{output_folder}{Specify the output folder for the summaries and the codebook. If no folder is specified, the function will return its output only to the R environment (unless this is suppressed under suppress_list_output).}

\item{suppress_list_output}{Logical. FALSE by default. If FALSE, the function will return a list output to the environment, containing the replication summaries and the codebook. If TRUE, these are not returned to the environment.}
}
\value{
The function create_replication_summaries returns a list consisting of two elements: A codebook and a list of data frames. Each data frame contains all replication summary statistics for the according replication (/effect).
The summary statistics returned (including their standard error) are the means and standard deviations for control and experimental groups, pooled standard deviations, raw mean differences and standardized mean differences (Hedge's g according to Borenstein et al., 2009).
}
\description{
\loadmathjax{}
\(
\\let\\underscore_
\)

Function to compute  replication aggregates from  person level data. Components of the standardized mean difference and their standard errors are calculated and reported. This is the first function (and the third computational step) of the MetaPipeX pipeline. For more details on the replication statistics, refer to the Details section. For more details on the pipeline, refer to the documentation of the MetaPipeX-package.
}
\details{
### Replication Statistics

All components of the standardized mean difference and their standard errors are returned by the function. Each standard error is returned to enable a meta-analysis on each component. The components and their standard errors are implemented as follows. Unless other sources are provided, effect size statistics are calculated according to Borenstein et al., 2009. The metafor::escalc function was used for SMD and MD (Viechtbauer, 2010). \cr
## mean (M)
\itemize{
 \item{R-Code} \cr
 \code{## apply the function} \cr
 \code{# treatment group mean (T_M):} \cr
 \code{mean(treatment_group$DV)} \cr
 \code{# control group mean (C_M):} \cr
 \code{mean(control_group$DV)} \cr
\item{Model} \cr
{
treatment group mean (T_M):
\mjdeqn{ \bar{x}_{T} = \frac{1}{n} \sum_{i \in T} x }{ \bar{x} \underscore{T} = \frac{1}{n} \sum \underscore{i \in T} x}
control group mean (C_M):
\mjdeqn{ \bar{x}_{C} = \frac{1}{n} \sum_{i \in C} x}{}
}
}
## standard error of the mean (SE_T_M, SE_C_M)
\itemize{
 \item{R-Code} \cr
 \code{## define the function} \cr
 \code{SE_of_mean_fct <- function(x)\{ } \cr
 \code{estimated_sd <- sqrt(sum((x-mean(x))^2)/(length(x)-1))} \cr
 \code{SE_of_mean <-  sd(x) / sqrt(length(x))} \cr
 \code{return(SE_of_mean)\} } \cr \cr
 \code{## apply the function} \cr
 \code{# standard error of treatment group mean (SE_T_M):} \cr
 \code{SE_of_mean_fct(treatment_group$DV)} \cr
 \code{# standard error of control group mean (SE_C_M):} \cr
 \code{SE_of_mean_fct(control_group$DV)} \cr
 \item{Model} \cr
 {
 \mjdeqn{ \hat{\sigma}_{\bar{x}} = \frac{\hat{\sigma}_{x}}{\sqrt{n}} = \sqrt{\frac{\frac{1}{n-1}\sum_{i=1}^n(x - \bar{x})^2}{n}} }{}
 }
}
## standard deviation (T_SD, C_SD)
\itemize{
\item{R-Code} \cr
\code{## apply the function} \cr
\code{# treatment group standard deviation (T_SD):} \cr
\code{ sd(treatment_group$DV)} \cr
\code{# control group standard deviation (C_SD):} \cr
\code{sd(control_group$DV)}
 \item{Model} \cr
{
\mjdeqn{ \hat{\sigma}  = \sqrt{ \frac{ \sum(x-\bar{x}^2) }{n-1}   } }{}
}
}
## standard error of the standard deviation (SE_T_SD, SE_C_SD)
\itemize{
\item{R-Code} \cr
\code{## define the function} \cr
\code{SE_SD_fct <- function(x)\{ } \cr
\code{SE_SD <- sd(x) / sqrt(2*(length(x)-1)) # for large n } \cr
\code{return(SE_SD) \} } \cr \cr
\code{## apply the function} \cr
\code{# standard error of the treatment group standard deviation (SE_T_SD):} \cr
\code{SE_SD_fct(treatment_group$DV)} \cr
\code{# standard error of the control group standard deviation (SE_C_SD):} \cr
\code{SE_SD_fct(control_group$DV)}
 \item{Model} \cr
{
\mjdeqn{ \hat{\sigma}_{\hat{\sigma}} = \frac{\hat{\sigma}_{x}}{\sqrt{2(n-1)}} = \sqrt{\frac{\frac{1}{n-1}\sum_{i=1}^n(x - \bar{x})^2}{2(n-1)}} }{}
\mjeqn{ \hat{\sigma}_{\hat{\sigma}} }{} is a simplified version of \mjeqn{ \sigma_{K_{n}S} }{} in Ahn & Fessler (2003). The authors demonstrate that for n > 10 it is reasonable to use Kn = 1. As for the overwhelming majority of samples n > k may be assumed, we excluded the term \mjeqn{K_{n}}{}.
}
}
## mean difference (MD)
\itemize{
\item{R-Code} \cr
\code{## apply the function} \cr
\code{metafor::escalc( } \cr
\code{measure = "MD", } \cr
\code{m1i = mean(treatment_group$DV),} \cr
\code{m2i = mean(control_group$DV), } \cr
\code{sd1i = sd(treatment_group$DV), } \cr
\code{sd2i = sd(control_group$DV), } \cr
\code{n1i = length(treatment_group$DV), } \cr
\code{n2i = length(control_group$DV), } \cr
\code{vtype = "HO" # assuming homoscedasticity } \cr
\code{)$yi } \cr
 \item{Model} \cr
{
\mjdeqn{ D = \bar{x}_{T} -  \bar{x}_{C} }{}
}
}
## standard error of mean difference (SE_MD)
\itemize{
\item{R-Code} \cr
\code{## apply the function} \cr
\code{metafor::escalc( } \cr
\code{measure = "MD", } \cr
\code{m1i = mean(treatment_group$DV),} \cr
\code{m2i = mean(control_group$DV), } \cr
\code{sd1i = sd(treatment_group$DV), } \cr
\code{sd2i = sd(control_group$DV), } \cr
\code{n1i = length(treatment_group$DV), } \cr
\code{n2i = length(control_group$DV), } \cr
\code{vtype = "HO" # assuming homoscedasticity } \cr
\code{)$vi } \cr
 \item{Model} \cr
{
\mjdeqn{ \hat{\sigma}_{\bar{x}_{T} -  \bar{x}_{C}} = \sqrt{ \frac{n_{T}+ n_{C}}{ n_{T} n_{C} } \sigma^2_{TC} } = \sqrt{ \frac{n_{T}+ n_{C}}{ n_{T} n_{C} } \frac{ \sum_{i = 1}^n (x_{T}-\bar{x}_{T})^2 + \sum_{i = 1}^n (x_{C}-\bar{x}_{C})^2 }{ n_{T} + n_{C} - 2 }   } }{}
}
}
## pooled standard deviation (pooled_SD)
\itemize{
\item{R-Code} \cr
\code{## define the function} \cr
\code{pooled_SD_fct <- function(t,c)\{ } \cr
\code{pooled_SD <- sqrt((} \cr
\code{(sum((t-mean(t))^2))+ # sample sum of squares sums treatment group} \cr
\code{(sum((c-mean(c))^2)) # sample sum of squares control group)/} \cr
\code{(length(t) + length(c) -2) # n+n-2} \cr
\code{) # end of sqrt} \cr
\code{return(pooled_SD)\} } \cr
\code{## apply the function} \cr
\code{pooled_SD_fct(treatment_group$DV, control_group$DV)}
\item{Model} \cr
{
\mjdeqn{ \hat{\sigma}_{TC} = \sqrt{ \frac{ \sum_{i = 1}^n (x_{T}-\bar{x}_{T})^2 + \sum_{i = 1}^n (x_{C}-\bar{x}_{C})^2 }{ n_{T} + n_{C} - 2 } } }{}
}
}
## standard error of pooled standard deviation (SE_pooled_SD)
\itemize{
\item{R-Code} \cr
\code{## define the function} \cr
\code{SE_pooled_SD_fct <- function(t,c)\{ } \cr
\code{pooled_SD <- sqrt((} \cr
\code{(sum((t-mean(t))^2))+ # sample sum of squares sums treatment group} \cr
\code{(sum((c-mean(c))^2)) # sample sum of squares control group)/} \cr
\code{(length(t) + length(c) -2) # n+n-2} \cr
\code{) # end of sqrt} \cr
\code{SE_pooled_SD <- pooled_SD/sqrt(2*(length(t)+length(c)-1))} \cr
\code{return(SE_pooled_SD)\} } \cr
\code{## apply the function} \cr
\code{SE_pooled_SD_fct(treatment_group$DV, control_group$DV)} \cr

\item{Model} \cr
{
\mjdeqn{ \hat{\sigma}_{\hat{\sigma}_{TC}} = \frac{ \hat{\sigma}_{TC} }{ \sqrt{ 2(n_{T}+n_{C}-1) } } }{}
The standard error is equivalent to that of the standard deviation. For further information, refer to the "standard error of the standard deviation" section.
}
}
## standardized mean difference (SMD)
\itemize{
\item{R-Code} \cr
\code{## apply the function} \cr
\code{metafor::escalc(} \cr
\code{measure = "SMD",} \cr
\code{m1i = mean(treatment_group$DV),} \cr
\code{m2i = mean(control_group$DV),} \cr
\code{sd1i = sd(treatment_group$DV),} \cr
\code{sd2i = sd(control_group$DV),} \cr
\code{n1i = length(treatment_group$DV),} \cr
\code{n2i = length(control_group$DV),} \cr
\code{vtype = "LS2" # Borenstein variance} \cr
\code{)$yi} \cr
\code{## apply the function} \cr
\item{Model} \cr
{
\mjdeqn{ g = d \left( 1 - \frac{3}{4(n_{T}+n_{C}-2) -1} \right)  }{}
with
\mjdeqn{ d =  \frac{ \bar{x}_{T} -  \bar{x}_{C} }{ \sqrt{ \frac{ \sum_{i = 1}^n (x_{T}-\bar{x}_{T})^2 + \sum_{i = 1}^n (x_{C}-\bar{x}_{C})^2 }{ n_{T} + n_{C} - 2 } }   }}{}
}
}
## standard error of standardized mean difference (SE_SMD)
\itemize{
\item{R-Code} \cr
\code{## apply the function} \cr
\code{sqrt(metafor::escalc(} \cr
\code{measure = "SMD",} \cr
\code{m1i = mean(treatment_group$DV),} \cr
\code{m2i = mean(control_group$DV),} \cr
\code{sd1i = sd(treatment_group$DV),} \cr
\code{sd2i = sd(control_group$DV),} \cr
\code{n1i = length(treatment_group$DV),} \cr
\code{n2i = length(control_group$DV),} \cr
\code{vtype = "LS2" # Borenstein variance} \cr
\code{)$vi)} \cr
\code{## apply the function} \cr
\item{Model} \cr
{
\mjdeqn{ \hat{\sigma}_{g} = \sqrt{ \hat{\sigma}_{d}^2 \left( 1 - \frac{3}{4(n_{T}+n_{C}-2) -1} \right)^2 } }{}
with
\mjdeqn{ \hat{\sigma}_{d}^2 = \frac{n_{T}+n_{C}}{n_{T}n_{C}} + \frac{d^2}{2(n_{T}+n_{C})} }{}
}
}
}
\examples{

# create vectors with names
# two multi-labs
MultiLab_names <- c("Multi_Lab_1", "Multi_Lab_2") # two projects
# two replication projects per multi lab
ReplicationProject_names <- c("Effect_A", "Effect_B", "Effect_C", "Effect_D")
# five replications per replication project
Replication_names <- c("Lab_A", "Lab_B", "Lab_C", "Lab_D", "Lab_E",
                       "Lab_A", "Lab_B", "Lab_C", "Lab_D", "Lab_E",
                       "Lab_F", "Lab_G", "Lab_H", "Lab_I", "Lab_J",
                       "Lab_F", "Lab_G", "Lab_H", "Lab_I", "Lab_J")

# create df with all example data
set.seed(1973)
example_data_df <- data.frame(MultiLab = rep(MultiLab_names, each = 100),
                              ReplicationProject = rep(ReplicationProject_names, each = 50),
                              # n = 10 (5 in control, 5 in treatment group)
                              Replication = rep(Replication_names, each = 10),
                              # random sampling for simulated data
                              DV = round(stats::rnorm(n = 2e2, mean = 0, sd = 5), 0),
                              Treatment = rep(c(1,0), times = 100))

# split the data per replication project to prepare for use in MetaPipeX::full_pipeline()
example_data_list <- split(example_data_df,
                           example_data_df$ReplicationProject)

## applying the input to the MetaPipeX function

# run create_replication_summaries
example_MetaPipeX_output <- MetaPipeXUpdate::create_replication_summaries(
data = example_data_list,
MultiLab = "MultiLab",
ReplicationProject = "ReplicationProject",
Replication = "Replication",
DV = "DV",
Group = "Treatment"
)

\dontrun{
All examples with additional comments are available on github:
https://github.com/JensFuenderich/MetaPipeX/tree/main/Supplementary_Material/Code_Examples
}

}
\references{
Ahn, S., & Fessler, J. A. (2003). Standard errors of mean, variance, and standard deviation estimators. EECS Department, The University of Michigan, 1(2).

Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R. (2009). Introduction to Meta-Analysis John Wiley & Sons. Ltd, Chichester, UK. 10.1002/9780470743386

Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. Journal of Statistical Software, 36(3), 1-48. doi: 10.18637/jss.v036.i03
}
