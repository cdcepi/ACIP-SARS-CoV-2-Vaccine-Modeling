#' @title Rescale incidence
#' @description Rescales incidence to match a given total number of infections. Can adjust to reduce incidence by a presumed background
#'    rate.
#' @param incidence Weekly incidence
#' @param totalInfections Total number of infections. It can represent the total number of infections or 
#'    the fraction of the population that became infected (totalInfections < 1).
#' @param weeksOfBackground Weeks of data at start and end of the data to presume reflect "background" incidence and to reduce the
#'    incidence curve by. Default to 3.
#' @return A vector of weekly incidence. If totalInfections represents a fraction of the population (< 1), then 
#'    this represents a rate of infection.
#' @keywords internal
#' @author Matt Clay <clay.matt@gmail.com>
#' @export
rescaleIncidence <- function(incidence, totalInfections, weeksOfBackground = 3) {

  if (sum(incidence < 0))
    stop("incidence cannot be < 0")

  if (totalInfections <= 0)
    stop("totalInfections cannot be <= 0")

  if ( (weeksOfBackground - round(weeksOfBackground) != 0) || weeksOfBackground < 0 || weeksOfBackground >= 0.5 * length(incidence))
    stop("weeksOfBackground must be an integer between 0 and 1/2 length(incidence)")

  # Take background to be the mean of the weeksOfBackground
  if (weeksOfBackground != 0)
    incidence.background <- mean(c(incidence[seq_len(weeksOfBackground)], 
                                   incidence[seq.int(length(incidence) - weeksOfBackground - 1, length(incidence))]))
  else
    incidence.background <- 0

  # Change to remove background rate
  incidence.rescaled <- pmax.int(incidence - incidence.background, 0)

  # Remove baseline incidence. Should likely make this more detailed by finding something more than min.
  # incidence.rescaled <- incidence - min(incidence)

  week <- seq_along(incidence)

  AUC <- sum(diff(week) * (incidence.rescaled[seq_len(length(incidence.rescaled) - 1)] + 
                           incidence.rescaled[seq.int(2, length(incidence.rescaled))]) / 2)
  incidence.rescaled <- incidence.rescaled / AUC

  # Return the incidence
  return(incidence.rescaled * totalInfections)
}
