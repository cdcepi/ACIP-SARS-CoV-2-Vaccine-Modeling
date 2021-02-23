#' @title Get vaccine uptake multiplier from overall vaccine allocation fractions
#' @description Calculates vaccine uptake multipliers suitable for use in flumodels::SEIRVModel. It takes
#'    vaccine allocation by age group (e.g., 30% to <18, 40% to 19-60, and 30% to 61+) and converts that into
#'    relative uptakes for arbitrary age bins.
#' @param ages Vector of ages. Each element represents the upper range of an age range. The lowest bound is presumed to be zero.
#'    If the oldest age does not reach the end of the population range, an additional element is added to span the full range.
#'    The final age bracket cannot start after the final bracket of originalContactMatrixAges.
#' @param population Population data frame that contains age and population values. Defaults to 2015 US population.
#'    The two required columns of the data frame are AGE and TOT_POP
#' @param year Year to sample population from. Defaults to 2015.
#' @param originalAllocation Original allocation fraction of vaccine, by age group.
#' @param originalAllocationAges Ages associated with originalAllocation. Must be in same format as ages.
#' @return Vector of vaccine uptake multipliers.
#' @author Matt Clay <clay.matt@gmail.com>
#' @importFrom dplyr filter select mutate summarise
#' @importFrom magrittr "%>%"
#' @export
getVaccineUptakeMultiplierFromAllocation <- function(ages,
                                                     population = flumodelsutil::flumodelsutil_data$population.US %>%
                                                       dplyr::filter(YEAR == year & MONTH == 1 & AGE != 999) %>%
                                                       dplyr::select(AGE, TOT_POP),
                                                     year = 2015,
                                                     originalAllocation, originalAllocationAges ) {

  #Check variables
  if(nrow(population) == 0)
    stop("Population not provided in proper format")

  if(sum(ages < 0))
    stop("ages cannot be < 0")
  if(sum(originalAllocationAges < 0))
    stop("originalAllocationAges cannot be < 0")
  if(sum(originalAllocation < 0))
    stop("originalAllocation cannot be < 0")

  if(is.unsorted(ages))
    stop("Ages must be increasing order")
  if (!all.equal(sum(originalAllocation), 1))
    stop("originalAllocation must sum to one")

  if (length(originalAllocationAges) != length(originalAllocation) &
      length(originalAllocationAges) != length(originalAllocation) -1 )
    stop("originalAllocationAges must be same length as originalAllocation")

  originalAllocationAgesFractions <- getPopulationFractions(originalAllocationAges)

  totalPop <- sum(population[,"TOT_POP"])

  # Add information on uptake to population data
  vaccineUsers <- population %>% dplyr::filter(AGE < originalAllocationAges[1]) %>%
                                  dplyr::mutate(useFraction = originalAllocation[1], popFraction = originalAllocationAgesFractions[1])

  if (length(originalAllocationAges) > 1) {
    for (age in 2:length(originalAllocationAges))
      vaccineUsers <- vaccineUsers %>% bind_rows (population %>% dplyr::filter(AGE < originalAllocationAges[age] &
                                                                          AGE >= originalAllocationAges[age-1]) %>%
                                                    dplyr::mutate(useFraction = originalAllocation[age],
                                                           popFraction = originalAllocationAgesFractions[age]))

    vaccineUsers <- vaccineUsers %>% bind_rows (population %>%
                                                  dplyr::filter(AGE >= originalAllocationAges[length(originalAllocationAges)]) %>%
                                                  dplyr::mutate(useFraction = originalAllocation[length(originalAllocation)],
                                                         popFraction = originalAllocationAgesFractions[length(originalAllocationAgesFractions)]) )
  }

  outputFractions <- numeric(length(ages) + 1)

  # Now compute relative uptake for new groups
  outputFractions[1] <- vaccineUsers %>% dplyr::filter(AGE <= ages[1]) %>%
    dplyr::mutate(use = TOT_POP * useFraction / popFraction) %>% dplyr::summarise(sum(use)) %>% unlist
  if (length(originalAllocationAges) > 1) {
    for (age in 2:length(ages))
      outputFractions[age] <- vaccineUsers %>% dplyr::filter(AGE <= ages[age] & AGE > ages[age-1]) %>%
        dplyr::mutate(use = TOT_POP * useFraction / popFraction) %>%
        dplyr::summarise(sum(use)) %>% unlist

    outputFractions[length(outputFractions)] <- vaccineUsers %>% dplyr::filter(AGE > ages[length(ages)]) %>%
      dplyr::mutate(use = TOT_POP * useFraction / popFraction) %>%
      dplyr::summarise(sum(use)) %>% unlist
  }

  # Rescale so it's relative on nice human-readable ratio
  outputFractions <- outputFractions / getPopulationFractions(ages) / totalPop


  return(outputFractions)

}
