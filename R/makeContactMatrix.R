#' @title Make contact matrix
#' @description Rebins an arbitrary contact matrix to a new group of population ranges. The default behavior is to use the UK POLYMOD data.
#' @param ages Vector of ages. Each element represents the upper range of an age range. The lowest bound is presumed to be zero.
#'    If the oldest age does not reach the end of the population range, an additional element is added to span the full range.
#'    The final age bracket cannot start after the final bracket of originalContactMatrixAges.
#' @param originalContactMatrix Contact matrix to serve as the basis for extrapolation. Defaults to UK POLYMOD.
#'    This must be a square matrix.
#' @param originalContactMatrixAges Age ranges associated with originalContactMatrix. Defaults to UK population in 5-year bins, with
#'    the final bracket starting at age 70. The format is a data frame with columns representing the age brackets and having two rows:
#'    AgeStart and AgeEnd. The first column has ageStart of 0. The last column has ageEnd of NA.
#' @param originalPopulationFractions Vector of age fractions associated with originalContactMatrixAges.
#'    Length must equal the dimension of originalContactMatrix. The vector will be normalized, so it can represent population fractions
#'    or population in each age bin.
#' @return A contact matrix.
#' @author Jason Asher <jason.m.asher@gmail.com>
#' @export
makeContactMatrix <- function( ages, originalContactMatrix = flumodels_data$POLYMOD.matrix, 
                               originalContactMatrixAges = flumodels_data$POLYMOD.age.ranges,
                               originalPopulationFractions = flumodels_data$population.fractions.US) {

  if (is.unsorted(ages))
    stop("Ages must be increasing order")

  if (sum(ages - round(ages)) != 0  | min(ages) < 1)
    stop("Ages must be positive integers")

  if (ages[length(ages)] > originalContactMatrixAges[1, length(originalContactMatrixAges)])
    stop(paste("Final age range is older than the maximum age possible based upon originalContactMatrixAges:",
               originalContactMatrixAges[1, length(originalContactMatrixAges)]))

  # Should check here to see if the matrix is one that we've determined before. This will help speed things up considerably.

  if(nrow(originalContactMatrix) != ncol(originalContactMatrix))
    stop("originalContactMatrix is not a square matrix")

  if (nrow(originalContactMatrix) != length(originalContactMatrixAges) ||
      length(originalContactMatrixAges) != length(originalPopulationFractions))
    stop("Dimension mismatch between originalContactMatrix, originalContactMatrixAges, and originalPopulationFractions")

  newAgeRanges <- makeAgeRangesFromInputs(ages)

  #Setting adjusts whether entered data is normalized by column (default) or row (assuming the original Mossong data has been transposed)
  byColumn <- TRUE

  #Calculate the matrix of expected number of contacts between single-year age groups
  if (byColumn) {
    matrixOfContacts <- t(apply(originalContactMatrix, 1, function(t){t*originalPopulationFractions})) #Multiply each column by the corresponding population
  } else {
    matrixOfContacts <- apply(originalContactMatrix, 2, function(t){t*originalPopulationFractions}) #Multiply each row by the corresponding population
  }

  #Symmetrize this matrix of contacts (to remove artifacts from study and to ensure consistency with the population)
  symmetrizedMatrixOfContacts <- (matrixOfContacts + t(matrixOfContacts))/2

  #Re-group the matrix of expected contacts
  regroupedMatrixOfContacts <- matrix(0, nrow = ncol(newAgeRanges), ncol = ncol(newAgeRanges),
                                      dimnames = list(names(newAgeRanges),names(newAgeRanges))) #Initialize with a zero matrix
  for (newRowIndex in seq_along(newAgeRanges)) {
    for (newColumnIndex in seq_along(newAgeRanges)){
      for (rowIndex in seq_along(originalContactMatrixAges)) {
        for (columnIndex in seq_along(originalContactMatrixAges))
          regroupedMatrixOfContacts[newRowIndex, newColumnIndex] <-
            regroupedMatrixOfContacts[newRowIndex, newColumnIndex] +
            (symmetrizedMatrixOfContacts[rowIndex, columnIndex] *
               getAgeRangeFraction(rowIndex, as.numeric(unlist(newAgeRanges[, newRowIndex])),
                                   ages = originalContactMatrixAges) *
               getAgeRangeFraction(columnIndex, as.numeric(unlist(newAgeRanges[, newColumnIndex])),
                                   ages = originalContactMatrixAges))
      }
    }
  }

  #debugOutput(regroupedMatrixOfContacts)

  #Re-normalize the new matrix of expected contacts
  newPopulation <- apply(newAgeRanges, 2, function(x){getPopulationForAgeRange(ageRange = as.numeric(unlist(x)),
                                                                               ages = originalContactMatrixAges,
                                                                               population = originalPopulationFractions)})

  if (byColumn) {
    newContactMatrix <- t(apply(regroupedMatrixOfContacts, 1, function(t){t/newPopulation})) #Divide each column by the corresponding population
  } else {
    newContactMatrix <- apply(regroupedMatrixOfContacts, 2, function(t){t/newPopulation}) #Divide each row by the corresponding population
  }

  return(newContactMatrix)
}



#Returns the population that would be assigned to the given age range c(ageStart, ageEnd)
# Ages is POLYMODAgeRanges
#Assumes constant interpolation between age groups when subdividing ranges
getPopulationForAgeRange <- function(ageRange, ages, population) {
  ageRangePopulation <- 0
  for (columnIndex in seq_along(ages)) {
    ageRangePopulation <- ageRangePopulation + (getAgeRangeFraction(columnIndex, ageRange, ages) * as.numeric(population[columnIndex]))
  }
  ageRangePopulation
}

#Returns the fraction of the corresponding age range with column index columnIndex that is between ageRange = c(ageStart,ageEnd)
# Ages is POLYMODAgeRanges
#Warning, cannot subdivide a column whose upper bound is NA - this function will return the whole population in that case
getAgeRangeFraction <- function(columnIndex, ageRange, ages) {
  ageStart <- ageRange[1]
  ageEnd <- ageRange[2]
  if (is.na(ages[2, columnIndex])) {
    if (is.na(ageEnd) || (ages[1, columnIndex] <= ageEnd)) {
      #Return
      1
    } else {
      #Return
      0
    }
  } else {
    overlapEnd <- min(ageEnd, ages[2, columnIndex], na.rm = TRUE)
    overlapStart <- max(ageStart, ages[1, columnIndex])
    yearsInOverlapOfAgeRanges <- max(overlapEnd - overlapStart + 1, 0)
    #Return
    yearsInOverlapOfAgeRanges / getAgeRangeLength(columnIndex, ages)
  }
}

getAgeRangeLength <- function(index, ages) {
  ages[2, index] - ages[1, index ] + 1
}


# Changes the age ranges provided into what the rebinning module prefers
makeAgeRangesFromInputs <- function(ages) {

  newFrame <- data.frame("col1" = c(AgeStart = 0, AgeEnd = ages[1]))
  names(newFrame) <- paste0("Age00to", ages[1])

  if (length(ages) != 1) {
    for (currentAge in seq.int(2, length(ages))) {
      newFrame[, paste0("Age", ages[currentAge-1]+1, "to", ages[currentAge])] <-
        c(AgeStart = ages[currentAge-1]+1, AgeEnd = ages[currentAge])
    }
  }
  newFrame[, paste0("Age", ages[length(ages)]+1, "plus")] <-
    c(AgeStart = ages[currentAge]+1, AgeEnd = NA)

  return(newFrame)

}
