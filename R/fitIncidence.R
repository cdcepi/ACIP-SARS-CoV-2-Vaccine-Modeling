#' @title Fit incidence
#' @description Fits a given incidence curve. Fits R0, generation time, seed infections, and prior immunity fractions (1 for each age group).
#'    Uses a model from the flumodels package. All parameters that are provided as a list are presumed
#'    to be something to fit. The list can be of either length two with the first list entry should be the low and the second list entry
#'    the high or of length three, where the first list entry is the low, the second is the start value for the fit,
#'    and the third is the high value. If no start value is provided, the start value for the fit is taken to be the average of the low
#'    and high values.
#'
#'    The only difference in parameters between this and the SEIRModel-family is that
#'    latentPeriod and infectiousPeriod are re-defined as generationTime and latentPercent. It is
#'    expected that a fit may fix latentPercent and then allow generationTime to vary rather than
#'    fitting both latent and infectious periods.
#'
#'    Only the required parameters from SEIRModel are listed below. Please see SEIRModel and the other models from flumodels
#'    for a full list of parameters.
#' @param incidence Data frame of weekly incidence. Columns are week & incidence, corresponding to week number [0-53] &
#'    incidence in population. Not applicable as a fit parameter. Incidence is overall, not broken down by population group.
#' @param population Size of population; defaults to 1
#' @param populationFractions Vector of population fractions (all non-negative, sum to 1); defaults to 1, 
#' representing a single population group
#' @param contactMatrix A matrix whose (row i, column j) entry denotes the number of potentially infectious 
#' contacts a single individual from group j has with individuals from group i each day; defaults to 
#' proportional mixing
#' @param R0 Average number of secondary cases from a single infected individual
#'   in a completely susceptible population; must be specified
#' @param generationTime Generation time in days; must be specified.
#' @param latentPercent Percent of total generation time spent in the latent period; must be specified.
#' @param seedInfections Number of infections to seed the outbreak with; must be specified and greater than 0.
#' @param method Optimiziation method. Options are L-BFGS-B from optim or any other from the list available in optim (which will cause
#'    the alabama package's constrOptim.nl method to be called). Not applicable as a fit parameter. Defaults to L-BFGS-B.
#' @param seedWeeksBeforeIncidence Weeks before the start of the supplied incidence curve to seed initial infections. Integer. Defaults to 0.
#' @param attackRatesToMatch Population-group specific attack rates to match. Not applicable as a fit parameter.
#' @param symptomaticIncidence Is the incidence serological or symptomatic. Defaults to false (serological).
#' @param fractionSymptomatic The fraction of cases that are symptomatic. Not a parameter to fit. Only used
#'     if symptomaticIncidence is set to true. Defaults to 0.5.
#' @param weightIncidenceTimeseries What relative weight to the penalty function should the incidence timeseries receive? Defaults to 1.
#' @param weightAttackRates What relative weight to the penalty function should the age-specific attack rates receive? Defaults to 1.
#' @param weightOverallAttackRate What relative weight to the penalty function should the overall attack rate receive? Defaults to 10.
#' @param ... Any other parameter to pass to a flumodel. Any of these may be specified as a list and thus fit. The parameters specified
#'    (vaccine, 2-dose vaccine, antiviral, etc.) will determine the type of flumodel used.
#' @return A fitIncidence object
#' @import flumodels
#' @importFrom alabama constrOptim.nl
#' @author Matt Clay <clay.matt@gmail.com>
#' @export
fitIncidence <- function(incidence, population, populationFractions, contactMatrix,
                         R0, generationTime, latentPercent, seedInfections = 1,
                         method = "L-BFGS-B", seedWeeksBeforeIncidence = 0, 
                         attackRatesToMatch, symptomaticIncidence = FALSE, fractionSymptomatic = 0.5, 
                         weightIncidenceTimeseries = 1, weightAttackRates = 1, weightOverallAttackRate = 10, ...) {
  
  if (missing(incidence))
    stop("Incidence must be specified")
  
  if(missing(attackRatesToMatch))
    stop("Attack rates to match must be specified")
  
  # Check parameters
  if (!is.data.frame(incidence))
    stop("Incidence must be a data frame")
  
  if (missing(generationTime))
    stop("Generation time must be specified")
  
  if (missing(latentPercent))
    stop("latentPercent must be specified")
  
  if (sum(seedInfections) <= 0)
    stop("seedInfections must be >0 to fit to anything.")
  
  # Make list of argument and ensure that it's in alphabetical order
  argument.list <- as.list(match.call())[-1]
  # argument.list <- argument.list[order(names(argument.list))]
  argument.list <- lapply(argument.list, eval)
  
  parameters.fixed <- argument.list[!unlist(lapply(argument.list, is.list))]
  
  # Skip these variables when thinking of what might be variable
  arguments.internal <- c("incidence", "method", "attackRatesToMatch", "seedWeeksBeforeIncidence", "fractionSymptomatic",
                          "weightIncidenceTimeseries", "weightAttackRates", "weightOverallAttackRate",
                          "symptomaticIncidence")
  parameters.SEIR.fixed <- parameters.fixed[ -which(names(parameters.fixed) %in% arguments.internal)]
  # Set the minimal length of the simulation to be 35 weeks or ( the amount of incidence data we have + the seed week offset )
  parameters.SEIR.fixed <- c(parameters.SEIR.fixed,
                             list(simulationLength = max(seedWeeksBeforeIncidence + nrow(incidence), 35) * 7 ))
  parameters.internal <- argument.list[which(names(argument.list) %in% arguments.internal)]
  argument.list <- argument.list[ -which(names(argument.list) %in% arguments.internal)]
  
  
  # Set parameters.internal to defaults if they weren't passed arguments
  if (sum(names(parameters.internal) == "symptomaticIncidence") == 0)
    parameters.internal$symptomaticIncidence <- symptomaticIncidence
  if (sum(names(parameters.internal) == "method") == 0)
    parameters.internal$method <- method
  if (sum(names(parameters.internal) == "fractionSymptomatic") == 0)
    parameters.internal$fractionSymptomatic <- fractionSymptomatic
  if (sum(names(parameters.internal) == "seedWeeksBeforeIncidence") == 0)
    parameters.internal$seedWeeksBeforeIncidence <- seedWeeksBeforeIncidence
  
  parameters.to.fit <- argument.list[unlist(lapply(argument.list, is.list))]
  
  if (length(parameters.to.fit) == 0)
    stop("Number of parameters to fit must not be zero")
  
  # Enforce that parameters to fit are of proper length
  if (sum(unlist(lapply(parameters.to.fit, length)) > 3) > 0 ||
      sum(unlist(lapply(parameters.to.fit, length)) < 2) > 0)
    stop("All parameters to fit must be lists of length two or three.")
  
  # Make upper and lower bounds
  bound.upper <- numeric()
  bound.lower <- numeric()
  bound.start <- numeric()
  # Deal with cases where start bounds were or were not specified
  # Use <<- notation to be able to affect outside of lapply
  appendParametersToBounds <- function(parameter) {
    bound.lower <<- c(bound.lower, unlist(parameter[[1]]))
    bound.upper <<- c(bound.upper, unlist(parameter[[length(parameter)]]))
    if (length(parameter) == 3)
      bound.start <<- c(bound.start, unlist(parameter[[2]]))
    else
      bound.start <<- c(bound.start, 0.5 * (parameter[[1]] + parameter[[2]]))
  }
  lapply(parameters.to.fit, appendParametersToBounds)
  
  # Check which parameters used so we know which SEIRModel to use
  argument.V.names <- c("vaccineAdministrationRatePerDay", "vaccineAvailabilityByDay", "vaccineUptakeMultiplier", "VEs", "VEp", "VEi",
                        "vaccineEfficacyDelay")
  argument.V2.names <- c("VEs1", "VEs2", "VEi1", "VEi2", "VEp1", "VEp2", "dose2Delay")
  
  argument.T.names <- c("fractionSeekCare", "fractionDiagnosedAndPrescribedOutpatient", "fractionAdhere",
                        "fractionAdmitted", "fractionDiagnosedAndPrescribedInpatient", "AVEi", "AVEp")
  if (sum(which(names(argument.list) %in% argument.T.names))) {
    # SEIRT-type models
    if (sum(which(names(argument.list) %in% argument.V2.names)))
      parameters.internal$SEIRModel.optimize <- SEIRModel.optimize <- "SEIRTV2DoseModel"
    else if (sum(which(names(argument.list) %in% argument.V.names)))
      parameters.internal$SEIRModel.optimize <- SEIRModel.optimize <- "SEIRTVModel"
    else
      parameters.internal$SEIRModel.optimize <- SEIRModel.optimize <- "SEIRTModel"
  } else if (sum(which(names(argument.list) %in% argument.V2.names))) {
    parameters.internal$SEIRModel.optimize <- SEIRModel.optimize <- "SEIRV2DoseModel"
  } else if (sum(which(names(argument.list) %in% argument.V.names))) {
    parameters.internal$SEIRModel.optimize <- SEIRModel.optimize <- "SEIRVModel"
  } else {
    parameters.internal$SEIRModel.optimize <- SEIRModel.optimize <- "SEIRModel"
  }
  
  # List of parameters to send to optim
  parameters.opt <- c(list(par = bound.start),
                      fn = optimize.SEIRV,
                      parameters.SEIR.fixed,
                      parameters.SEIR.fixed.names = list(names(parameters.SEIR.fixed)),
                      parameters.to.fit,
                      parameters.to.fit.names = list(names(parameters.to.fit)),
                      parameters.internal,
                      list(upper = bound.upper),
                      list(lower = bound.lower),
                      list(control = list(factr = 1e7))
  )
  
  # parscale seems quite off, leading the optimizer to totally incorrect results.
  #, parscale = bound.start
  
  if (method == "L-BFGS-B")
    opt <- do.call("optim", parameters.opt)
  
  else {
    constraint.jacobian <- constraint.jacobian.create(variables = length(bound.fit))
    # parameters.opt$control.optim <- parameters.opt$control
    names(parameters.opt[names(parameters.opt) == "control"]) <- "control.optim"
    opt <- do.call("alabama::constrOptim.nl", c(parameters.opt, list(hin = constraint.function,
                                                            hin.jac = constraint.jacobian.function,
                                                            constraint.jacobian = constraint.jacobian)))
  }
  
  bound.fit <- opt$par
  
  # Now to translate bound.fit into parameter list for SEIRModel
  parameters.fitted <- vector(mode = "list", length = length(parameters.to.fit))
  names(parameters.fitted) <- names(parameters.to.fit)
  fitLocation <- 1
  for (parameterIterator in seq_along(parameters.to.fit)) {
    parameter <- names(parameters.to.fit)[parameterIterator]
    value <- parameters.to.fit[[which(names(parameters.to.fit) == parameter)]]
    value.length <- length(value[[1]])
    parameters.fitted[[parameterIterator]] <- unlist(bound.fit[fitLocation:(fitLocation + value.length - 1)])
    fitLocation <- fitLocation + value.length
  }
  
  # Hack to deal with generation time
  if (sum(names(parameters.fitted) == "generationTime") > 0 & sum(names(parameters.fitted) == "latentPercent") > 0) {
    # Both were fitted
    parameters.fitted$latentPeriod <- parameters.fitted$generationTime * parameters.fitted$latentPercent
    parameters.fitted$infectiousPeriod <- parameters.fitted$generationTime * (1 - parameters.fitted$latentPercent)
    parameters.fitted$generationTime <- NULL
    parameters.fitted$latentPercent <- NULL
  } else if (sum(names(parameters.fitted) == "generationTime") > 0) {
    # Only generation time was fit
    parameters.fitted$latentPeriod <- parameters.fitted$generationTime * parameters.fixed$latentPercent
    parameters.fitted$infectiousPeriod <- parameters.fitted$generationTime * (1 - parameters.fixed$latentPercent)
    parameters.fitted$generationTime <- NULL
    parameters.SEIR.fixed$latentPercent <- NULL
  } else if (sum(names(parameters.fitted) == "latentPercent") > 0) {
    # Only latentPercent was fit
    parameters.fitted$latentPeriod <- parameters.fixed$generationTime * parameters.fitted$latentPercent
    parameters.fitted$infectiousPeriod <- parameters.fixed$generationTime * (1 - parameters.fitted$latentPercent)
    parameters.fitted$latentPercent <- NULL
    parameters.SEIR.fixed$generationTime <- NULL
  } else {
    # Both fixed
    parameters.fitted$latentPeriod <- argument.list$generationTime * argument.list$latentPercent
    parameters.fitted$infectiousPeriod <- argument.list$generationTime * (1 - argument.list$latentPercent)
    parameters.SEIR.fixed$latentPercent <- NULL
    parameters.SEIR.fixed$generationTime <- NULL
  }
  
  parameters.total <- c(parameters.fitted,
                        parameters.SEIR.fixed)
  
  # Manually remove latentPercent
  parameters.total$latentPercent <- NULL
  
  flumodels.wrapper <- function(fun) {
    get(fun, asNamespace("flumodels"))
  }
  
  model <- do.call(flumodels.wrapper(SEIRModel.optimize), parameters.total)
  
  data <- list(model = model,
               fit = opt,
               parameters.fitted = parameters.fitted,
               parameters.fixed = parameters.fixed,
               parameters.total = parameters.total,
               parameters.to.fit = parameters.to.fit,
               parameters.internal = parameters.internal)
  
  class(data) <- "fitIncidence"
  return(data)
}


#' @title Compute fit value for SEIR models
#' @description Computes the fit value
#' @param bound.fit Guess for current iteration.
#' @param incidence Overall incidence curve to match
#' @param attackRatesToMatch Final attack rates (by age group) to match
#' @param SEIRModel.optimize SEIRModel function that we should use to simulate outbreak.
#' @param parameters.to.fit.names Names of parameters to fit.
#' @param parameters.SEIR.fixed.names Names of fixed parameters for SEIR model.
#' @param symptomaticIncidence Is the incidence serological or symptomatic. Defaults to false (serological).
#' @param fractionSymptomatic The fraction of cases that are symptomatic. Not a parameter to fit. Only used
#'     if symptomaticIncidence is set to true. Defaults to 0.5.
#' @param seedWeeksBeforeIncidence Days before the start of the supplied incidence curve to seed initial infections.
#' @param weightIncidenceTimeseries What relative weight to the penalty function should the incidence timeseries receive? Defaults to 1.
#' @param weightAttackRates What relative weight to the penalty function should the age-specific attack rates receive? Defaults to 1.
#' @param weightOverallAttackRate What relative weight to the penalty function should the overall attack rate receive? Defaults to 10.
#' @return Fit value for optim
#' @keywords internal
optimize.SEIRV <- function(bound.fit, incidence, attackRatesToMatch, SEIRModel.optimize, 
                           parameters.to.fit.names,
                           parameters.SEIR.fixed.names, 
                           symptomaticIncidence = FALSE,
                           fractionSymptomatic = 0.5,
                           seedWeeksBeforeIncidence, 
                           weightIncidenceTimeseries = 1,
                           weightAttackRates = 1,
                           weightOverallAttackRate = 10, ...) {
  
  argument.list <- as.list(match.call())[-1]
  
  parameters.to.model <- vector(mode = "list", length = length(parameters.to.fit.names))
  names(parameters.to.model) <- parameters.to.fit.names
  fitLocation <- 1
  for (parameterIterator in seq_along(parameters.to.fit.names)) {
    value <- argument.list[[which(names(argument.list) == parameters.to.fit.names[parameterIterator])]]
    value.length <- length(value[[1]])
    parameters.to.model[[parameterIterator]] <- unlist(bound.fit[fitLocation:(fitLocation + value.length - 1)])
    fitLocation <- fitLocation + value.length
  }
  # Hack to deal with generation time
  if (sum(names(parameters.to.model) == "generationTime") > 0 & sum(names(parameters.to.model) == "latentPercent") > 0) {
    # Both were fitted
    parameters.to.model$latentPeriod <- parameters.to.model$generationTime * parameters.to.model$latentPercent
    parameters.to.model$infectiousPeriod <- parameters.to.model$generationTime * (1 - parameters.to.model$latentPercent)
    parameters.to.model$generationTime <- NULL
    parameters.to.model$latentPercent <- NULL
  } else if (sum(names(parameters.to.model) == "generationTime") > 0) {
    # Only generation time was fit
    parameters.to.model$latentPeriod <- parameters.to.model$generationTime * argument.list$latentPercent
    parameters.to.model$infectiousPeriod <- parameters.to.model$generationTime * (1 - argument.list$latentPercent)
    parameters.to.model$generationTime <- NULL
    parameters.SEIR.fixed.names <- parameters.SEIR.fixed.names[-which(parameters.SEIR.fixed.names == "latentPercent")]
  } else if (sum(names(parameters.to.model) == "latentPercent") > 0) {
    # Only latentPercent was fit
    parameters.to.model$latentPeriod <- argument.list$generationTime * parameters.to.model$latentPercent
    parameters.to.model$infectiousPeriod <- argument.list$generationTime * (1 - parameters.to.model$latentPercent)
    parameters.to.model$latentPercent <- NULL
    parameters.SEIR.fixed.names <- parameters.SEIR.fixed.names[-which(parameters.SEIR.fixed.names == "generationTime")]
  } else {
    # Both fixed
    parameters.to.model$latentPeriod <- argument.list$generationTime * argument.list$latentPercent
    parameters.to.model$infectiousPeriod <- argument.list$generationTime * (1 - argument.list$latentPercent)
    parameters.SEIR.fixed.names <- parameters.SEIR.fixed.names[-which(parameters.SEIR.fixed.names == "generationTime")]
    parameters.SEIR.fixed.names <- parameters.SEIR.fixed.names[-which(parameters.SEIR.fixed.names == "latentPercent")]
  }
  
  incidence <- incidence$incidence
  
  parameters.SEIR.fixed <- argument.list[which(names(argument.list) %in% parameters.SEIR.fixed.names)]
  
  parameters.to.model <- c(parameters.to.model,
                           parameters.SEIR.fixed)
  
  model <- do.call(SEIRModel.optimize, parameters.to.model)
  
  if (!("population" %in% names(parameters.to.model)))
    parameters.to.model$population <- 1
  
  SEIR.penalty.function(model = model, incidence = incidence, 
                        attackRatesToMatch = attackRatesToMatch, populationFractions = parameters.to.model$populationFractions,
                        population = parameters.to.model$population,
                        symptomaticIncidence = symptomaticIncidence,
                        fractionSymptomatic = fractionSymptomatic,
                        seedWeeksBeforeIncidence = seedWeeksBeforeIncidence,
                        weightIncidenceTimeseries = weightIncidenceTimeseries,
                        weightAttackRates = weightAttackRates,
                        weightOverallAttackRate = weightOverallAttackRate)
}

#' @title Penalty function for SEIR-type models
#' @description Computes the penalty compared to preexisting incidence.
#' @param model SEIR-type model
#' @param incidence Array of weekly incidence
#' @param attackRatesToMatch Attack rates for population groups
#' @param populationFractions Population fractions
#' @param population Total population, Default: 1
#' @param symptomaticIncidence Is the incidence serological or symptomatic. Defaults to false (serological).
#' @param fractionSymptomatic The fraction of cases that are symptomatic. Not a parameter to fit. Only used
#'     if symptomaticIncidence is set to true. Defaults to 0.5.
#' @param seedWeeksBeforeIncidence Days before the start of the supplied incidence curve to seed initial infections.
#' @param weightIncidenceTimeseries What relative weight to the penalty function should the incidence timeseries receive? Defaults to 1.
#' @param weightAttackRates What relative weight to the penalty function should the age-specific attack rates receive? Defaults to 1.
#' @param weightOverallAttackRate What relative weight to the penalty function should the overall attack rate receive? Defaults to 10.
#' @return Total penalty
#' @keywords internal
#' @importFrom flumodels getInfectionTimeSeries getInfections
SEIR.penalty.function <- function(model, incidence, attackRatesToMatch, populationFractions, 
                                  population = 1, symptomaticIncidence = FALSE, fractionSymptomatic = 0.5,
                                  seedWeeksBeforeIncidence = 0, weightIncidenceTimeseries = 1, 
                                  weightAttackRates = 1, weightOverallAttackRate = 30) {
  infections.by.week.overall <- flumodels::getInfectionTimeSeries(model, incidence = TRUE, byGroup = FALSE, asRate = TRUE,
                                                                  symptomatic = symptomaticIncidence, fractionSymptomatic = fractionSymptomatic, 
                                                                  byWeek = TRUE)
  infections.by.week.overall <- infections.by.week.overall[is.na(infections.by.week.overall) == FALSE]
  
  # Simple L2 norm of differences in overall infections (rescaled to be a rate)
  infections.by.week.error <- sqrt(sum( ( infections.by.week.overall[seq_along(incidence) + seedWeeksBeforeIncidence] -
                                         incidence / population )^2 ))  / max(incidence / population)
  
  # Then compare attack rates by age category
  infections <- flumodels::getInfections(model, asRate = TRUE, symptomatic = symptomaticIncidence, fractionSymptomatic = fractionSymptomatic)
  infections.error <- sqrt(sum( ((infections / populationFractions - attackRatesToMatch) / max(attackRatesToMatch))^2))
  
  # Now compare overall attack rate
  infections.overall.error <- abs(sum(infections) - sum(attackRatesToMatch * populationFractions)) / max(attackRatesToMatch)
  
  return ( (weightIncidenceTimeseries * infections.by.week.error  + weightAttackRates * infections.error + 
              weightOverallAttackRate * infections.overall.error) / 
             (weightIncidenceTimeseries + weightAttackRates + weightOverallAttackRate))
}


#' @title Constraint function
#' @description Computes the constraint function for use in constrOptim.nl (hin parameter)
#' @return Constraint function for upper / lower bounds
#' @keywords internal
constraint.function <- function(x, lower, upper, ...) {
  constraint <- numeric(2*length(x))
  
  for (currentVariable in seq_along(x)) {
    constraint[(currentVariable-1)*2+1] <- x[currentVariable] - lower[currentVariable]
    constraint[currentVariable*2] <- upper[currentVariable] - x[currentVariable]
  }
  
  return(constraint)
  
}


#' @title Creates the Jacobian of the constraint function
#' @description Creates Jacobian of the constraint function for use in constrOptim.nl
#'    (hin.jac parameter)
#' @return Jacobian of constraint function
#' @keywords internal
constraint.jacobian.create <- function(x, variables = length(x)) {
  jacobian <- matrix(0, 2 * variables, variables)
  
  for (currentVariable in seq_along(variables)) {
    jacobian[(currentVariable-1)*2+1,] <- c( rep(0, currentVariable - 1), 1,
                                             rep(0, variables - currentVariable) )
    jacobian[currentVariable*2,] <- c( rep(0, currentVariable - 1), -1,
                                       rep(0, variables - currentVariable) )
  }
  
  return(jacobian)
}


#' @title Jacobian of the constraint function
#' @description Returns the previously-computed Jacobian of the constraint function for
#'    use in constrOptim.nl (hin.jac parameter)
#' @return Jacobian of constraint function
#' @keywords internal
constraint.jacobian.function <- function(constraint.jacobian, ...) {
  return (constraint.jacobian)
}
