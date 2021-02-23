#' @name flumodelsutil_data
#' @title Internal and external data for flumodelsutil
#' @description A collection of data for use with \code{flumodelsutil}. Most are used by internal functions but are exposed here for use by
#' external parties. Those items most likely to be used externally are \code{pandemic.deaths} and \code{population.US.total}.
#' @docType data
#' @format A \code{list} of data. The named entries are described below:
#' \describe{
#'\item{ILI.data}{
#'A \code{data.frame} with the following columns:
#' \itemize{
#' \item REGION TYPE: Only value is National
#' \item REGION: Only value is NA (since it's national data)
#' \item YEAR: Integers with values 1997 to 2015 corresponding to calendar year
#' \item Season: Integers with values 1997 to 2015 corresponding to the start of the year X/X+1 flu season.
#' \item WEEK: ILINet week
#' \item Date: Start date of the week.
#' \item \% WEIGHTED ILI: ILINet reported % of weighted hospital cases due to ILI.
#' \item \%UNWEIGHTED ILI: ILINet reported % of unweighted hospital cases due to ILI.
#' \item AGE 0-4: Number of ILI cases in the 0-4 age bracket.
#' \item AGE 25-49: Number of ILI cases in the 25-49 age bracket.
#' \item AGE 25-64: Number of ILI cases in the 24-64 age bracket.
#' \item AGE 5-24: Number of ILI cases in the 5-24 age bracket.
#' \item AGE 50-64: Number of ILI cases in the 50-64 age bracket.
#' \item AGE 65: Number of ILI cases in the 64+ age bracket.
#' \item ILITOTAL: Total number of ILI reported cases.
#' \item NUM. OF PROVIDERS: Number of providers reporting to ILINet.
#' \item TOTAL PATIENTS: Total number of patients reported to ILINet.
#' }}
#' 
#' \item{pandemic.deaths}{
#' A \code{data.frame} with 63 rows and 5 variables:
#' \describe{
#'   \item{\code{year}}{Year corresponding to data. \code{integer}}
#'   \item{\code{week}}{CDC ILI week. \code{integer}}
#'   \item{\code{week.of.pandemic}}{How many weeks since the beginning of the pandemic. Starts at 1 for the first recorded element of the given pandemic. \code{integer}}
#'   \item{\code{pandemic}}{Year of pandemic. Included pandemics are 1918, 1957, 1968, 2009. \code{integer}}
#'   \item{\code{deaths}}{Estimated mortality. Please note that not all estimates are for the entire US. \code{integer}
#'   \describe{
#'   \item{1918 data}{Excess death rates (annual basis) per 100k from influenza and pneumonia in 47 US cities.}
#'   \item{1957 data}{Influenza and pneumonia deaths in the United States.}
#'   \item{1968 data}{Influenza and pneumonia  deaths in 122 cities in the United States.}
#'   \item{2009 data}{Influenza deaths in the United States. This includes both waves. For wave 2 only, just use 2009 week 31 (week.of.pandemic 18) and onwards.}
#'   }}
#'}}
#'
#'\item{POLYMOD.age.ranges}{
#'A \code{data.frame}, where each column corresponds to an age range. The
#' column names denote the age ranges. The two rows are AgeStart and AgeEnd.
#' The given age range corresponds to individuals in the range [AgeStart, AgeEnd].
#' A value of NA for AgeEnd corresponds to no upper bound, so [70,NA] corresponds
#' to all individuals aged 70 or more.}
#' 
#' \item{POLYMOD.matrix}{
#' A \code{double} matrix, matrix where k_ij is proportional to the observed number of contacts
#' (both physical and nonphysical) that a respondent in age band j makes with other
#' individuals in age band i. The names of both columns and rows denote the
#' age ranges.}
#' 
#' \item{population.UK}{
#' An \code{integer} array, one entry per age range.
#' The \code{names} denote the age ranges.}
#' 
#' \item{population.US}{
#' A \code{data.frame} with many columns. Those of interest for this package are:
#' \itemize{
#' \item YEAR: Which year the data corresponds to. Currently only 2015 data is provided.
#' \item AGE: Age of the population of interest. Age of X corresponds to [X,X+1),
#' so 0 means infants under one year of age. Age of 999 corresponds to the entire population.
#' \item TOT_POP: Total population in the given age bin (or entire US population if AGE == 999).}}
#' 
#' \item{population.US.total}{
#' A \code{numeric} of total 2015 US population.}
#' }
#' 
#' @source  
#' \describe{
#' \item{ILI.data}{
#' CDC ILINet}
#' \item{pandemic.deaths}{
#' Mossong J, Hens N, Jit M, Beutels P, Auranen K, et al. (2008)
#' Social Contacts and Mixing Patterns Relevant to the Spread of Infectious Diseases.
#' PLOS Medicine 5(3): e74. https://doi.org/10.1371/journal.pmed.0050074}
#' 
#' \item{POLYMOD.age.ranges}{
#' Mossong J, Hens N, Jit M, Beutels P, Auranen K, et al. (2008)
#' Social Contacts and Mixing Patterns Relevant to the Spread of Infectious Diseases.
#' PLOS Medicine 5(3): e74. https://doi.org/10.1371/journal.pmed.0050074}
#' 
#' \item{POLYMOD.matrix}{
#' Mossong J, Hens N, Jit M, Beutels P, Auranen K, et al. (2008)
#' Social Contacts and Mixing Patterns Relevant to the Spread of Infectious Diseases.
#' PLOS Medicine 5(3): e74. https://doi.org/10.1371/journal.pmed.0050074}
#' 
#' \item{population.US}{
#' US Census}
#' 
#' \item{population.US.total}{
#' US Census}
#' }
"flumodelsutil_data"
