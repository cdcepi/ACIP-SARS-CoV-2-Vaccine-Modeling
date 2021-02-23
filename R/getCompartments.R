#' @title Get compartments
#' @description This is a utility function to obtain the names of the compartments with the specified type
#' @param model The model for which compartment names are desired
#' @param type The type of the compartment
#' @return A vector of compartment names
#' @keywords internal
getCompartments <- function(model, type) {
  UseMethod("getCompartments", model)
}

#' @rdname getCompartments
#' @method getCompartments SEIRModel
#' @keywords internal
#' @export
getCompartments.SEIRModel <- function(model, type) {
  return(switch(type,
                S = "S",
                E = "E",
                A = "A",
                I = "I",
                R = "R"))
}

#' @rdname getCompartments
#' @method getCompartments SEIRVModel
#' @keywords internal
#' @export
getCompartments.SEIRVModel <- function(model, type) {
  return(switch(type,
                S          = c("S", "Sv"),
                E          = c("E", "Ev"),
                A          = c("A", "Av"),
                I          = c("I", "Iv"),
                R          = c("R", "Rv"),
                vaccinated = c("V")))
}

#' @rdname getCompartments
#' @method getCompartments SEAIRTVModel
#' @keywords internal
#' @export
getCompartments.SEAIRTVModel <- function(model, type) {
  return(switch(type,
                S          = c("S", "Sv","Svb"),
                E          = c("E", "Ev","Evb"),
                A          = c("A", "Av","Avb"),
                I          = c("I", "Iv","Ivb"),
                R          = c("R", "Rv","Rvb"),
                vaccinated = c("V")))
}

#' @rdname getCompartments
#' @method getCompartments SEIRV2DoseModel
#' @keywords internal
#' @export
getCompartments.SEIRV2DoseModel <- function(model, type) {
  return(switch(type,
                S          = c("S", "Sv", "Svb"),
                E          = c("E", "Ev", "Evb"),
                A          = c("A", "Av", "Avb"),
                I          = c("I", "Iv", "Ivb"),
                R          = c("R", "Rv", "Rvb"),
                vaccinated = c("V"))) #Only tracks 1st dose
}

#' @rdname getCompartments
#' @method getCompartments SEIRVMonoModel
#' @keywords internal
#' @export
getCompartments.SEIRVMonoModel <- function(model, type) {
  return(switch(type,
                S          = c("S", "Sv", "SvM"),
                E          = c("E", "Ev", "EvM"),
                I          = c("I", "Iv", "IvM"),
                R          = c("R", "Rv", "RvM"),
                vaccinated = c("V"))) #Only tracks first dose (which it tries to give to everyone)
}

#' @rdname getCompartments
#' @method getCompartments SEIRVPrimeBoostModel
#' @keywords internal
#' @export
getCompartments.SEIRVPrimeBoostModel <- function(model, type) {
  return(getCompartments.SEIRV2DoseModel(model, type))
}

#' @rdname getCompartments
#' @method getCompartments SEAIR_NPI_Model
#' @keywords internal
#' @export
getCompartments.SEAIR_NPI_Model <- function(model, type) {
  return(switch(type,
                S = c("S"),
                E = c("E.SymptomaticInfector","E.AsymptomaticInfector"),
                A = c("A.preQuarantine.SymptomaticInfector","A.preQuarantine.AsymptomaticInfector","A.Quarantine","A.nonQuarantine"),
                I = c("I.preIsolation", "I.Isolation","I.nonIsolation", "I.Quarantine.Symptomatic","I.Quarantine.Asymptomatic", "I.postQuarantine.Asymptomatic"),
                R = c("R")))
 
}

