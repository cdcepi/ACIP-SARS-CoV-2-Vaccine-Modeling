#' @title Get incidence from ILI data
#' @description Get incidence from CDC-collected ILI data
#' @param year Vector of ages. Each element represents the upper range of an age range. The lowest bound is presumed to be zero.
#'    If the oldest age does not reach the end of the population range, an additional element is added to span the full range.
#' @param ILI Data frame containing ILI data. Defaults to requesting data directly from ILINet via cdcfluview package.
#'    Required columns are year, week, and weighted_ili
#' @param startWeekOfSeason Week of the year (integer) that is considered the start of the flu season. Defaults to 40.
#' @param attackRateToMatch The overall attack rate for the season. Incidence will be rescaled to match this value.
#' @param population Population to scale the resulting incidence curve to. Defaults to 1.
#' @param weeksOfBackground Weeks of data at start and end of the data to presume reflect "background" incidence and to reduce the
#'    incidence curve by. Default to 3.
#' @return A data frame of weeks and incidence. It will only return weeks for which data exists and so may not correspond to a full year.
#'    Columns are: week & incidence.
#' @importFrom dplyr filter select rename
#' @importFrom cdcfluview get_flu_data
#' @importFrom magrittr "%>%"
#' @author Matt Clay <clay.matt@gmail.com>
#' @export
# getIncidenceFromILI <- function(year, ILI = ILI.data, startWeekOfSeason = 40, attackRateToMatch, population = 1, weeksOfBackground = 3) {
getIncidenceFromILI <- function(year, ILI = cdcfluview::ilinet(region = "national", years = (year-1):(year+1)),
                                startWeekOfSeason = 40, attackRateToMatch, population = 1, weeksOfBackground = 3) {

  if(nrow(ILI) == 0)
    stop("No ILI data")

  if((year - round(year)) != 0)
    stop("Year must be an integer")

  if((startWeekOfSeason - round(startWeekOfSeason)) != 0 | startWeekOfSeason < 0 | startWeekOfSeason > 53)
    stop("Year must be an integer on [0,53]")

  if (population <= 0)
    stop("Population must be > 0")

  if((weeksOfBackground - round(weeksOfBackground) != 0) || weeksOfBackground < 0 || weeksOfBackground >= 0.5*nrow(ILI))
    stop("weeksOfBackground must be an integer between 0 and 1/2 nrow(ILI)")

  data.year <- ILI %>% 
    dplyr::rename(ili.year = year) %>%
    dplyr::filter(( ili.year == year & week >= startWeekOfSeason) | (ili.year == (year + 1) & week < startWeekOfSeason)) %>%
    dplyr::select(week, `weighted_ili`) %>% dplyr::rename(ILI = `weighted_ili`) %>% dplyr::filter(!is.na(ILI))

  if (nrow(data.year) == 0)
    stop("No data for year specified")

  return ( data.frame(week = data.year$week,
             incidence = rescaleIncidence(data.year$ILI, attackRateToMatch * population,
                                          weeksOfBackground = weeksOfBackground) ) )
}
