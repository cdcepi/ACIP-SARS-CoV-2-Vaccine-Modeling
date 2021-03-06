% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SEIRVPrimeBoostModel.R
\name{checkInputs.SEIRVPrimeBoost}
\alias{checkInputs.SEIRVPrimeBoost}
\title{Check SEIR+V prime-boost inputs}
\usage{
checkInputs.SEIRVPrimeBoost(
  population,
  populationFractions,
  contactMatrix,
  R0,
  latentPeriod,
  infectiousPeriod,
  seedInfections,
  priorImmunity,
  useCommunityMitigation,
  communityMitigationStartDay,
  communityMitigationDuration,
  communityMitigationMultiplier,
  vaccineAdministrationRatePerDay,
  vaccineAvailabilityByDayPrime,
  vaccineAvailabilityByDayBoost,
  vaccineUptakeMultiplier,
  boostDelay,
  VEs1,
  VEs2,
  VEi1,
  VEi2,
  VEp1,
  VEp2,
  vaccineEfficacyDelay,
  simulationLength,
  seedStartDay,
  tolerance,
  method
)
}
\value{
List of parameters for the SEIR+V prime-boost model
}
\description{
Checks the input parameters for the SEIR+V prime-boost model
}
\keyword{internal}
