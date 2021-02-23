#' @title Get vaccine uptake multiplier from overall vaccine coverage fractions
#' @description Calculates vaccine uptake multipliers suitable for use in flumodels::SEIRVModel. It takes
#'    vaccine coverage by age group (e.g., 65% for <18, 15% for 19-60, and 55% for 61+) and converts that into
#'    relative uptakes for arbitrary age bins.
#' @param ages Vector of ages. Each element represents the upper range of an age range. The lowest bound is presumed to be zero.
#'    If the oldest age does not reach the end of the population range, an additional element is added to span the full range.
#'    The final age bracket cannot start after the final bracket of originalContactMatrixAges.
#' @param population Population data frame that contains age and population values. Defaults to 2015 US population.
#'    The two required columns of the data frame are AGE and TOT_POP
#' @param year Year to sample population from. Defaults to 2015.
#' @param originalCoverage Original coverage fraction of vaccine, by age group.
#' @param originalCoverageAges Ages associated with originalCoverage Must be in same format as ages.
#' @return Vector of vaccine uptake multipliers.
#' @author Matt Clay <clay.matt@gmail.com>
#' @importFrom dplyr filter select mutate summarise
#' @importFrom magrittr "%>%"
#' @export
getVaccineUptakeMultiplierFromCoverage <- function(ages,
                                                     population = flumodelsutil::flumodelsutil_data$population.US %>%
                                                       dplyr::filter(YEAR == year & MONTH == 1 & AGE != 999) %>%
                                                       dplyr::select(AGE, TOT_POP),
                                                     year = 2015,
                                                     originalCoverage, originalCoverageAges ) {

  #Check variables
  if(nrow(population) == 0)
    stop("Population not provided in proper format")

  if(sum(ages < 0))
    stop("ages cannot be < 0")
  if(sum(originalCoverageAges < 0))
    stop("originalCoverageAges cannot be < 0")
  if(sum(originalCoverage < 0) | sum(originalCoverage > 1))
    stop("originalCoverage cannot be < 0 or > 1")

  if(is.unsorted(ages))
    stop("Ages must be increasing order")

  if (length(originalCoverageAges) != length(originalCoverage) &
      length(originalCoverageAges) != length(originalCoverage) -1 )
    stop("originalCoverageAges must be same length as originalCoverage")

  originalCoverageAgesFractions <- getPopulationFractions(originalCoverageAges)

  # Compute the original allocation of vaccine and then call getVaccineUptakeMultiplierFromAllocation
  totalVaccinesUsedByAges <- originalCoverageAgesFractions * originalCoverage
  
  vaccineAllocation <- totalVaccinesUsedByAges / sum(totalVaccinesUsedByAges)

  return(getVaccineUptakeMultiplierFromAllocation(ages = ages,
                                                  population = population,
                                                  year = year,
                                                  originalAllocation = vaccineAllocation,
                                                  originalAllocationAges = originalCoverageAges))
}
