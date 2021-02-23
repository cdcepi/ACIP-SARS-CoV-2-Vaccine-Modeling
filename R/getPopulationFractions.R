#' @title Get population fractions
#' @description Gets population fractions for a given set of age ranges
#' @param ages Vector of ages. Each element represents the upper range of an age range. The lowest bound is presumed to be zero.
#'    If the oldest age does not reach the end of the population range, an additional element is added to span the full range.
#'    This function by default uses \code{flumodelsutil_data$population.US} for age information (2015 US population).
#'    One could also use the \code{acs} package to download data dynamically, but that requires the user to request an individual
#'    key from the ACS for access.
#' @param year Year to sample population from. Defaults to 2015.
#' @param population Population data frame that contains age and population values. Defaults to 2015 US population.
#'    The two required columns of the data frame are AGE and TOT_POP.
#' @return A vector population fractions that sums to 1
#' @author Matt Clay <clay.matt@gmail.com>
#' @export
getPopulationFractions <- function(ages,
                                   year = 2020,
                                   population = flumodels_data$population.US) {
  
  if (sum(names(match.call()) == "population") == 0){ 
    # population <- flumodelsutil::flumodelsutil_data$population.US[flumodelsutil::flumodelsutil_data$population.US$AGE != 999 & 
    #                                                               flumodelsutil::flumodelsutil_data$population.US$YEAR == year & 
    #                                                               flumodelsutil::flumodelsutil_data$population.US$MONTH == 1, c("AGE", "TOT_POP")]
  
      
      # population <- flumodels_data$population.US[flumodels_data$population.US$AGE != 999 & 
      #                                                               flumodels_data$population.US$YEAR == year & 
      #                                                               flumodels_data$population.US$MONTH == 1, c("AGE", "TOT_POP")]
}  
  
  population = 
    filter(population, 
           AGE!=999, 
           YEAR == year,
           MONTH == 1) %>% 
    select("AGE","TOT_POP")
  
  
  if(nrow(population) == 0)
    stop("No population for given year")

  if(is.unsorted(ages))
    stop("Ages must be increasing order")

  if(sum(ages - round(ages)) != 0  | min(ages) < 1)
    stop("Ages must be positive integers")

  maxAge <- max(population[,"AGE"])

  if (max(ages) >= maxAge)
    stop("Final age bracket includes no people. The final element in ages is the starting point of the final age bin. There is no need to include an endpoint representing the max age in the population.")

  ages <- append(ages, maxAge + 1)

  population.fractions <- vapply(ages, function(x) {sum(population[population$AGE <= x, "TOT_POP"])},
                                 FUN.VALUE = numeric(1))
  population.fractions <- c(population.fractions[1],
                            population.fractions[2:length(population.fractions)] - population.fractions[1:(length(population.fractions)-1)] ) /
    sum(population[,"TOT_POP"])
  names(population.fractions) <- c(paste0("0-", ages[1]),
    if(length(ages) > 1) {
           vapply(2:length(ages), function(x) { paste0(ages[x-1]+1, ifelse(ages[x] >= maxAge, "+", paste0("-", ages[x])))}, FUN.VALUE = character(1))
    }
  )

  return(population.fractions)
}
