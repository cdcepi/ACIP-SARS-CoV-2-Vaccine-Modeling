#' @title SEAIR+TV Model
#' @description A model for influenza that uses a SEAIR framework with age 
#'   structure and that allows for antiviral treatment and vaccination
#' @param population Size of population; defaults to 1
#' @param populationFractions Vector of population fractions (all non-negative, 
#'   sum to 1); defaults to 1, representing a single population group
#' @param contactMatrix A matrix whose (row i, column j) entry denotes the 
#'   number of potentially infectious contacts a single individual from group j 
#'   has with individuals from group i each day; defaults to proportional mixing
#' @param contactMatrices A list of contact matrices which can be specified instead of 
#' contactMatrix, and uses the same formatting. If contactMatrices is specified,
#' the parameter ContactMatrixStartDay must also be specified, to indicate the time period that
#' each contact matrix corresponds to.  
#' @param ContactMatrixStartDay a vector of length equal to contactMatrices. 
#' The first entry should be 0. each entry of ContactMatrixStartDay corresponds to the first day that
#' the corresponding entry of contactMatrices goes into effect. Increasing integer values only. 
#' @param R0 Average number of secondary cases from a single infected individual
#'   in a completely susceptible population; must be specified
#' @param latentPeriod Latent period in days; must be specified
#' @param infectiousPeriod Infectious period in days; must be specified
#' @param fractionLatentThatIsInfectious Fraction of latent period that is infectious; 
#'   in the range 0 to 1 (inclusive), must be specified
#' @param relativeInfectivityAsymptomatic Infectivity of asymptomatic 
#'   individuals, relative to symptomatic individuals; defaults to 1
#' @param seedInfections Fraction of the population to seed with infections; 
#'   single fraction or vector of fractions by population group; defaults to 0
#' @param priorImmunity Fraction of the population with prior immunity; single 
#'   fraction, or vector of fractions by population group; defaults to 0
#' @param useCommunityMitigation Whether or not to use community mitigation
#'   implemented by modulation of the contact matrix; defaults to FALSE
#' @param useSecondCommunityMitigation Whether or not to use a second community mitigation period
#'   implemented by modulation of the contact matrix; defaults to FALSE. Should only be used if useCommunityMitigation is true.
#' @param communityMitigationStartDay If using community mitigation, day of the
#'   simulation on which to start mitigation; must be specified if applicable
#' @param communityMitigationDuration If using community mitigation, duration of
#'   time during which mitigation is in effect; must be specified if applicable
#' @param secondCommunityMitigationDuration If using second community mitigation, duration of
#'   time during which second mitigation is in effect; must be specified if applicable
#' @param communityMitigationMultiplier If using community mitigation, the
#'   non-negative matrix of multipliers that will be used to modulate the contact
#'   matrix by elementwise multiplication; must be specified if applicable
#' @param secondCommunityMitigationMultiplier If using second community mitigation, the
#'   non-negative matrix of multipliers that will be used to modulate the contact
#'   matrix by elementwise multiplication; must be specified if applicable
#' @param fractionSymptomatic Fraction of the infections that are symptomatic 
#'   cases; single fraction or vector of fractions by population group;  defaults
#'   to 0.5
#' @param fractionSeekCare Fraction of the symptomatic cases that seek medical
#'   care; single fraction or vector of fractions by population group; defaults
#'   to 0.6
#' @param fractionDiagnosedAndPrescribedOutpatient Fraction of the outpatient
#'   medical care seeking cases that are diagnosed and prescribed antiviral
#'   drugs; single fraction or vector of fractions by population group; defaults
#'   to 0.7
#' @param fractionAdhere Fraction of the cases that are prescribed antiviral
#'   drugs that adhere to the regimen; single fraction or vector of fractions by
#'   population group; defaults to 0.8
#' @param fractionAdmitted Fraction of the cases that require hospitalization
#'   that are admitted; single fraction or vector of fractions by population
#'   group; defaults to 1
#' @param fractionDiagnosedAndPrescribedInpatient Fraction of the hospitalized
#'   cases that are diagnosed and prescribed antiviral drugs; single fraction or
#'   vector of fractions by population group; defaults to 1
#' @param AVEi Antiviral efficacy: prevention of transmission from infected 
#'   individuals taking antiviral drugs; defaults to 0
#' @param AVEp Antiviral efficacy: probability that antiviral treatment averts 
#'   hospitalization and/or death; defaults to 0
#' @param vaccineAdministrationRatePerDay Vaccine administration rate each day; 
#'   defaults to 0
#' @param vaccineAvailabilityByDayPrime Vector that contains the amount of prime vaccine 
#'   available each day; defaults to 0
#' @param vaccineAvailabilityByDayBoost Vector that contains the amount of boost vaccine 
#'   available each day; defaults to 0
#' @param vaccineUptakeMultiplier Vector of multipliers that determines the
#'   relative rate at which vaccine is given to each age group; defaults to
#'   vaccine being allotted proportionally by population
#' @param VEs1 Vaccine efficacy (first dose): prevention of infection for vaccinated susceptible 
#'   individuals; single fraction or vector of fractions by population group;
#'   defaults to 0
#' @param VEi1 Vaccine efficacy (first dose): prevention of transmission from vaccinated 
#'   infected individuals; single fraction or vector of fractions by population
#'   group; defaults to 0
#' @param VEp1 Vaccine efficacy (first dose): prevention of symptomatic illness in
#'   vaccinated infected individuals; single fraction or vector of fractions by population
#'   group; defaults to 0
#' @param VEs2 Vaccine efficacy (second dose): prevention of infection for vaccinated susceptible 
#'   individuals; single fraction or vector of fractions by population group;
#'   defaults to 0
#' @param VEi2 Vaccine efficacy (second dose): prevention of transmission from vaccinated 
#'   infected individuals; single fraction or vector of fractions by population
#'   group; defaults to 0
#' @param VEp2 Vaccine efficacy (second dose): prevention of symptomatic illness in
#'   vaccinated infected individuals; single fraction or vector of fractions by population
#'   group; defaults to 0
#' @param vaccineEfficacyDelay Delay in days between administration of dose and 
#'   onset of protection; defaults to 7
#' @param boostDelay Minimum days between prime and boost doses; defaults to 14
#' @param simulationLength Number of days to simulate after seeding infections; 
#'   defaults to 240
#' @param seedStartDay Day on which to seed initial infections; defaults to 0
#' @param tolerance Absolute tolerance for numerical integration; defaults to 
#'   1e-8
#' @param hospDuration a numeric vector of the mean duration of hospitalization, in days, for each age group.
#' @param caseHospitalizationRatio a numeric vector of the fraction of symptomatic cases in each population group that will require hospitalization.
#' @param DelayOnsetToHosp a single numeric value for the mean delay from symptom onset to hospitalization.
#' @param Seasonality a logical value for whether or not to implement seasonality in transmissibility. Defaults to false. If set to true, 
#' R0 will vary sinusoidally over a 365 day period, with the timing and strength of seasonal variation determined by the parameters dayOfR0,
#' offset, and R0_amp. 
#' @param offset the day at which seasonally varying R0 reaches its maximum value, relative to the start of the modeled period. 
#' Can be positive or negative (for example, if the first day of the modeled period is treated as January 1st, and transmission is assummed to peak
#' on December 31st, offset=-1 and offset=365 will produce the same seasonal pattern). Set to -3.8*7 default, indicating transmission efficiency 
#' peaking in early December, if the model starts on January 1st.  
#' @param R0_amp the difference in the maximum and minimum values of the seasonally varying R0. Set to 0 by default.
#' @param dayOfR0 the day, relative to the start of the modeled period, at which the seasonally varying R0 is equal to the given value
#' of the R0 parameter
#' @param vaccineUptakeMultipliers a more flexible alternative to the vaccineUptakeMultiplier, allowing the first dose in a vaccination series to be targeted
#' to different population groups during different time periods (second doses are distributed between population groups on a 
#' per-capita basis according to the number of eligible members (received the first dose at least boostDays ago)).
#' A list of vectors of the same format as the vaccineUptakeMultiplier, described above. If vaccineUptakeMultipliers is specified,
#' the parameter vaccineUptakeMultipliers must also be specified, to indicate the time period that
#' each vector of uptake multipliers corresponds to.  
#' @param vaccineUptakeMultipliers a vector of length equal to vaccineUptakeMultipliers 
#' The first entry should be 0. each entry of ContactMatrixStartDay corresponds to the first day that
#' the corresponding entry of vaccineUptakeMultipliers goes into effect. Increasing integer values only. 
#' @param method Which integration method to use. Defaults to lsoda
#' @return a SEAIRTVModel object
#' @export
SEAIRTVModel <- function(population, populationFractions, contactMatrix, R0,
                        latentPeriod, infectiousPeriod, seedInfections, priorImmunity,
                        fractionLatentThatIsInfectious, relativeInfectivityAsymptomatic,
                        useCommunityMitigation, communityMitigationStartDay,
                        communityMitigationDuration, communityMitigationMultiplier,
                        fractionSymptomatic,  fractionSeekCare, fractionDiagnosedAndPrescribedOutpatient,
                        fractionAdhere, fractionAdmitted, fractionDiagnosedAndPrescribedInpatient, AVEi, AVEp,
                        vaccineAdministrationRatePerDay, vaccineAvailabilityByDayPrime,
                        vaccineAvailabilityByDayBoost,
                        vaccineUptakeMultiplier, VEs1, VEs2, VEi1, VEi2, VEp1, VEp2, vaccineEfficacyDelay,
                        boostDelay=28, simulationLength, seedStartDay, tolerance, method,
                        useSecondCommunityMitigation=F, secondCommunityMitigationDuration, 
                        secondCommunityMitigationMultiplier, hospDuration, caseHospitalizationRatio, DelayOnsetToHosp,
                        Seasonality=F,offset=-3.8*7,R0_amp=0, dayOfR0=75, 
                        ContactMatrices, ContactMatrixStartDay=0,vaccineUptakeMultipliers=list(1), 
                        vaccineUptakeMultiplierStartDay=c(0)) {
  #Check inputs #TODO: Add checks for all inputs
  specifiedArguments <- names(match.call())[-1]

  argumentList <- lapply(specifiedArguments, as.name)
  names(argumentList) <- specifiedArguments
  if(!"contactMatrix" %in% names(argumentList)){argumentList$contactMatrix=ContactMatrices[[1]]}
  if(!"vaccineUptakeMultiplier" %in% names(argumentList) )
    {argumentList$vaccineUptakeMultiplier=vaccineUptakeMultipliers[[1]]}
  if(!"vaccineUptakeMultipliers" %in% names(argumentList) )
    {argumentList$vaccineUptakeMultipliers=vaccineUptakeMultipliers}
  if(!"vaccineUptakeMultiplierStartDay" %in% names(argumentList) )
    {argumentList$vaccineUptakeMultiplierStartDay=vaccineUptakeMultiplierStartDay}
  if(!"boostDelay" %in% names(argumentList) )
  {argumentList$boostDelay=boostDelay}
  
  parameters <- do.call("checkInputs.SEAIRTV", argumentList) #Get parameters from checked inputs
  
  parameters$useSecondCommunityMitigation = useSecondCommunityMitigation
  parameters$Seasonality = Seasonality
  
  if(! "boostDelay" %in% names(parameters)){stop("boostDelay not in names(parameters)")}
  if(Seasonality)
  {
    if(missing(offset)|missing(dayOfR0)|missing(R0_amp))
    {stop("offset, R0_amp, and dayOfR0 must be specified when transmission seasonality is used")}
  else
    {
    parameters$offset = offset
    parameters$dayOfR0 = dayOfR0
    parameters$R0_amp = abs(R0_amp)
    }  }
  if (useSecondCommunityMitigation) {
    if (missing(secondCommunityMitigationDuration)) {
      stop("secondCommunityMitigationDuration must be specified when using second community mitigation.", 
           call. = FALSE)
    }
    else
      {parameters$secondCommunityMitigationDuration = secondCommunityMitigationDuration  }
    if (missing(secondCommunityMitigationMultiplier)) {
      stop("secondCommunityMitigationMultiplier must be specified when using second community mitigation.", 
           call. = FALSE)
    }
    else
      {checkNonNegative(secondCommunityMitigationMultiplier)
        parameters$secondCommunityMitigationMultiplier = secondCommunityMitigationMultiplier}
    if (!all(dim(secondCommunityMitigationMultiplier) == dim(contactMatrix))) {
      stop("Dimensions of secondCommunityMitigationMultiplier do not match those of contactMatrix", 
           call. = FALSE)
    }
    
  }
  

  initialState <- with(parameters, {
    c(S  = (1 - priorImmunity) * populationFractions,
      E  = 0 * populationFractions,
      A  = 0 * populationFractions,
      I  = 0 * populationFractions,
      R  = priorImmunity * populationFractions,
      Sv = 0 * populationFractions,
      Ev = 0 * populationFractions,
      Av = 0 * populationFractions,
      Iv = 0 * populationFractions,
      Rv = 0 * populationFractions,
      Svb = 0 * populationFractions,
      Evb = 0 * populationFractions,
      Avb = 0 * populationFractions,
      Ivb = 0 * populationFractions,
      Rvb = 0 * populationFractions,
      V  = 0 * populationFractions,
      Vb = 0 * populationFractions,
      vaccinatingPrime = rep(1, length(populationFractions)),
      vaccinatingBoost = rep(0, length(populationFractions)),
      Symp = 0 * populationFractions,
      Hosp = 0 * populationFractions,
      Discharged = 0 * populationFractions)
  })
  
  rootFunction <- function(t, state, parameters) {
    stateList <- reconstructState.SEAIRV(state)
    with(append(stateList, parameters), {
      return(c(ifelse(vaccinatingPrime > 0, populationFractions - V - tolerance, 1),
               ifelse(vaccinatingBoost > 0, V - Vb - tolerance, 1),
               ifelse(V - Vb - tolerance > tolerance & vaccinatingBoost <= 0, 0, 1)))
    })
  }
  
  eventFunction <- function(t, state, parameters) {
    stateList <- reconstructState.SEAIRV(state)
    with(append(stateList, parameters), {
      state[getLabels("vaccinatingPrime", length(populationFractions))] <-
        ifelse(populationFractions - V > tolerance, 1, 0)
      state[getLabels("vaccinatingBoost", length(populationFractions))] <-
        ifelse(V - Vb > tolerance, 1, 0)
      return(state)
    })
  }
  
  
  rawOutput <- integrateModel(initialState = initialState,
                              parameters = parameters,
                              derivativeFunction = getDerivative.SEAIRTV,
                              seedFunction = doSeed.SEAIRV,
                              rootFunction = rootFunction,
                              eventFunction = eventFunction)
  
  #Build the SEAIRVModel object to return
  model <- list(parameters = parameters,
                rawOutput = rawOutput)
  class(model) <- c("SEAIRTVModel")
  return(model)
}

#' @title Check SEAIR+TV inputs
#' @description Checks the input parameters for the SEAIR+V model
#' @return List of parameters for the SEAIR+V model
#' @keywords internal
checkInputs.SEAIRTV <- function(population, populationFractions, contactMatrix, R0,
                               latentPeriod, infectiousPeriod, seedInfections, priorImmunity,
                               fractionLatentThatIsInfectious, relativeInfectivityAsymptomatic,
                               useCommunityMitigation, communityMitigationStartDay,
                               communityMitigationDuration, communityMitigationMultiplier,
                               fractionSymptomatic,  fractionSeekCare, fractionDiagnosedAndPrescribedOutpatient,
                               fractionAdhere, fractionAdmitted, fractionDiagnosedAndPrescribedInpatient, AVEi, AVEp,
                               vaccineAdministrationRatePerDay, vaccineAvailabilityByDayPrime,
                               vaccineAvailabilityByDayBoost,
                               vaccineUptakeMultiplier, VEs1, VEs2, VEi1, VEi2, VEp1, VEp2, vaccineEfficacyDelay,
                               boostDelay=28,
                               simulationLength, seedStartDay, tolerance, method, 
                               useSecondCommunityMitigation, secondCommunityMitigationDuration, 
                               secondCommunityMitigationMultiplier, hospDuration, caseHospitalizationRatio,
                               DelayOnsetToHosp, Seasonality=F,offset=-3.8*7,R0_amp=0,dayOfR0=75,
                               ContactMatrices, ContactMatrixStartDay,
                               vaccineUptakeMultipliers,vaccineUptakeMultiplierStartDay) 
  {
  specifiedArguments <- names(match.call())[-1]
  argumentList <- lapply(specifiedArguments, as.name)
  names(argumentList) <- specifiedArguments
  if(missing(contactMatrix) & (missing(ContactMatrices) | missing(ContactMatrixStartDay)))
  {stop("if contactMatrix isn't provided, both ContactMatrices and ContactMatrixStartDay must be")}
  if(missing(contactMatrix)){argumentList$contactMatrix = ContactMatrices[[1]]}

  if(!missing(vaccineAdministrationRatePerDay) & missing(vaccineUptakeMultiplier) & (missing(vaccineUptakeMultipliers) | missing(vaccineUptakeMultiplierStartDay)))
  {stop("if vaccinating in model, and vaccineUptakeMultiplier isn't provided, both vaccineUptakeMultipliers and vaccineUptakeMultiplierStartDay must be")}
  if(missing(vaccineUptakeMultiplier) & !missing(vaccineAdministrationRatePerDay)){argumentList$vaccineUptakeMultiplier = vaccineUptakeMultipliers[[1]]}
  SEIRParameters <- do.call("checkInputs.SEIR", argumentList)

  
  #fractionLatentThatIsInfectious
  if (missing(fractionLatentThatIsInfectious)) {
    stop("fractionLatentThatIsInfectious must be specified.", call. = FALSE)
  }
  checkBetween0and1(fractionLatentThatIsInfectious)
  #relativeInfectivityAsymptomatic
  checkPositive(relativeInfectivityAsymptomatic)
  
  if(!(missing(ContactMatrices) | missing(ContactMatrixStartDay)))
  {
    SEIRParameters$ContactMatrices = ContactMatrices
    SEIRParameters$ContactMatrixStartDay = ContactMatrixStartDay
  }
  
  
  

    for(k in 1:length(vaccineUptakeMultipliers))
      {
      checkNonNegative(vaccineUptakeMultipliers[[k]])
      checkDimensionsMatch(vaccineUptakeMultipliers[[k]], populationFractions)
      }
    SEIRParameters$vaccineUptakeMultipliers = vaccineUptakeMultipliers
    SEIRParameters$vaccineUptakeMultiplierStartDay = vaccineUptakeMultiplierStartDay
    #SEIRParameters$vaccinationRateAgeMultipliers = list()
    
    # for(k in 1:length(vaccineUptakeMultipliers))
    #   {
    #   SEIRParameters$vaccinationRateAgeMultipliers[[k]] = vaccineUptakeMultipliers[[k]]*populationFractions
    #   if(sum(SEIRParameters$vaccinationRateAgeMultipliers[[k]])>0)
    #     {
    #     SEIRParameters$vaccinationRateAgeMultipliers[[k]] = SEIRParameters$vaccinationRateAgeMultipliers[[k]]/sum(SEIRParameters$vaccinationRateAgeMultipliers[[k]])
    #     }
    #   else {
    #     warning("vaccineUptakeMultipliers prevents vaccination from occurring.", call. = FALSE)
    #   }
    #   }
  

  
  antiviralParameters <- do.call("checkInputs.Antiviral", argumentList)
  #Update arguments passed to checkInputs.Vaccine using SEIRParameters
  argumentList$population <- SEIRParameters$population
  argumentList$populationFractions <- SEIRParameters$populationFractions
  argumentList$seedStartDay <- SEIRParameters$seedStartDay
  argumentList$simulationLength <- SEIRParameters$simulationLength
  #vaccineParameters <- do.call("checkInputs.Vaccine", argumentList)
  vaccineParameters <- do.call("checkInputs.VaccinePrimeBoost", argumentList)
  #vaccineParameters <- do.call("checkInputs.Vaccine2Dose", argumentList)
  
  # Update SEIRParameters for beta, lambda1, and lambda2
  SEIRParameters$lambda1 = 1 / ((1-fractionLatentThatIsInfectious) * latentPeriod)
  SEIRParameters$lambda2 = 1 / (fractionLatentThatIsInfectious * latentPeriod)
  SEIRParameters$fractionLatentThatIsInfectious = fractionLatentThatIsInfectious
  SEIRParameters$relativeInfectivityAsymptomatic = relativeInfectivityAsymptomatic
  SEIRParameters$hospDuration = hospDuration
  SEIRParameters$caseHospitalizationRatio = caseHospitalizationRatio
  SEIRParameters$DelayOnsetToHosp = DelayOnsetToHosp
  SEIRParameters$boostDelay = boostDelay
  # If there are no community mitigations, eg contact matrix among presymptomatic class is the same as among symptomatic class
  SEIRParameters$beta = R0 / 
    max(Mod(eigen(
        (infectiousPeriod + 1 / SEIRParameters$lambda2) * contactMatrix,
      symmetric = FALSE, only.values = TRUE
    )$values))
  
  SEIRParameters$beta_amp = R0_amp / 
    max(Mod(eigen(
      (infectiousPeriod + 1 / SEIRParameters$lambda2) * contactMatrix,
      symmetric = FALSE, only.values = TRUE
    )$values))
  #Return the parameters
  
  if(! "boostDelay" %in% names(SEIRParameters)){stop("boostDelay not in names(vaccineParameters)")}
  return(c(SEIRParameters, antiviralParameters, vaccineParameters))
}

#This function implements the multivariate derivative of the SEAIR+TV model
#parameters should define populationFractions, contactMatrix, beta, lambda, gamma,
#AVEi.eff, VEs, VEi, and the function vaccinationRate(t)
#Note that the total population is normalized to be 1
getDerivative.SEAIRTV <- function(t, state, parameters) {
  stateList <- reconstructState.SEAIRV(state)
  with(append(stateList, parameters), {
    
    #print(t)
    
    if(Seasonality)
    {
      day = t + 1
      if(length(day)!=1){stop(paste0("length day != 1, length day = ",length(day)))}
      if(length(beta)!=1){stop(paste0("length beta != 1, length beta = ",length(beta)))}
      if(length(R0_amp)!=1){stop(paste0("length R0_amp != 1, length R0_amp = ",length(R0_amp)))}
      if(length(offset)!=1){stop(paste0("length offset != 1, length offset = ",length(offset)))}
      
      beta_Max =  beta + beta_amp*0.5 - 0.5*beta_amp*cos((2*pi/(52*7))*(dayOfR0-(offset)))
      beta_time = 0.5*beta_amp*cos((2*pi/(52*7))*(day-(offset))) + (beta_amp/2)+beta_Max-beta_amp
      if(length(beta_time)!=1){stop(paste0("length beta_time != 1, length beta_time = ",length(beta_time)))}
      beta = beta_time
    }
    
    
    
    if (useCommunityMitigation) {
      if ((t >= communityMitigationStartDay) && (t < communityMitigationEndDay)) {
        contactMatrix <- communityMitigationMultiplier * contactMatrix
      } 
    }
    if(useSecondCommunityMitigation)
    {
      if((t>communityMitigationEndDay) & (t <= communityMitigationEndDay+secondCommunityMitigationDuration))
      {
        contactMatrix <- secondCommunityMitigationMultiplier * contactMatrix
      }
    }
    
      contactMatrix = ContactMatrices[[max(which(ContactMatrixStartDay<=t))]]
     
      # replacing negative counts with 0
      S = (abs(S)+S)/2 
      Sv = (abs(Sv)+Sv)/2 
      Svb = (abs(Svb)+Svb)/2 
      
      E = (abs(E)+E)/2 
      Ev = (abs(Ev)+Ev)/2 
      Evb = (abs(Evb)+Evb)/2 
      
      A = (abs(A)+A)/2 
      Av = (abs(Av)+Av)/2 
      Avb = (abs(Avb)+Avb)/2 
      
      I = (abs(I)+I)/2 
      Iv = (abs(Iv)+Iv)/2 
      Ivb = (abs(Ivb)+Ivb)/2 
      
      R = (abs(R)+R)/2 
      Rv = (abs(Rv)+Rv)/2 
      Rvb = (abs(Rvb)+Rvb)/2 
      
      V = (abs(V)+V)/2 
      Vb = (abs(Vb)+Vb)/2 
      

      EligibleForDose1 = ifelse(populationFractions - V < 0, 0 , populationFractions - V)
      EligibleForDose2 = ifelse(V - Vb < 0, 0 , V - Vb) 

      #if(any(EligibleForDose1<0)){stop(paste("any(EligibleForDose1<0), t =",t,", group =",which(EligibleForDose2<0)))}
      #if(any(EligibleForDose2<0)){stop(paste("any(EligibleForDose2<0), t =",t,", group =",which(EligibleForDose2<0)))}
      
      vaxMultiplier = vaccineUptakeMultipliers[[max(which(vaccineUptakeMultiplierStartDay<=max(t- vaccineEfficacyDelay + 1,0)))]]*populationFractions
      if(sum(vaxMultiplier)>0){vaxMultiplier = vaxMultiplier/sum(vaxMultiplier)}
      else{warning("vaccineUptakeMultipliers prevents vaccination from occurring.", call. = FALSE)}
      if(any(vaxMultiplier<0)){stop(paste("any(vaxMultiplier<0), t =",t))}
      
      isVaccinatingPrimeByAge <- (vaccinatingPrime > 0) & (populationFractions - V > 0)
      effectiveVaccinationPrimeMultiplier <- sum(ifelse(isVaccinatingPrimeByAge, 1, 0) * vaxMultiplier)
      if (effectiveVaccinationPrimeMultiplier > 0) {
        vaccinationRatePrimeByAge <- vaccinationRatePrime(t) * vaxMultiplier / effectiveVaccinationPrimeMultiplier
      } else {vaccinationRatePrimeByAge <- 0}
      

      vaxMultiplierBoost = rep(1,length(populationFractions))*EligibleForDose2
      if(sum(vaxMultiplierBoost)>0){vaxMultiplierBoost = vaxMultiplierBoost/sum(vaxMultiplierBoost)}
      vaxMultiplierBoost = 
        vaccineUptakeMultipliers[[max(which(vaccineUptakeMultiplierStartDay<=max(t- vaccineEfficacyDelay - boostDelay + 1,0)))]]*vaxMultiplierBoost
      if(sum(vaxMultiplierBoost)>0){vaxMultiplierBoost = vaxMultiplierBoost/sum(vaxMultiplierBoost)}
      isVaccinatingBoostByAge <- (vaccinatingBoost > 0) & (V - Vb > 0)
      effectiveVaccinationBoostMultiplier <- sum(ifelse(isVaccinatingBoostByAge, 1, 0) * vaxMultiplierBoost )
      if (effectiveVaccinationBoostMultiplier > 0) {
        vaccinationRateBoostByAge <- vaccinationRateBoost(t ) * vaxMultiplierBoost /   effectiveVaccinationBoostMultiplier
      } else {vaccinationRateBoostByAge <- 0}

  
    
      if(any(vaccinationRateBoostByAge<0)){stop(any(vaccinationRateBoostByAge<0))}
     
    

    # effectiveVaccinationDose1Multiplier <- sum(ifelse(V < populationFractions, 1, 0) * vaccinationRateAgeMultiplier)
    # if (effectiveVaccinationDose1Multiplier > 0) {
    #   vaccinationRateDose1ByAge <- vaccinationRateDose1(t) * vaccinationRateAgeMultiplier /
    #     effectiveVaccinationDose1Multiplier
    # } else {
    #   vaccinationRateDose1ByAge <- 0
    # }
    # 
    # effectiveVaccinationDose2Multiplier <- sum(ifelse(Vb < V, 1, 0) * vaccinationRateAgeMultiplier)
    # if (effectiveVaccinationDose2Multiplier > 0) {
    #   vaccinationRateDose2ByAge <- vaccinationRateDose2(t) * vaccinationRateAgeMultiplier /
    #     effectiveVaccinationDose2Multiplier
    # } else {
    #   vaccinationRateDose2ByAge <- 0
    # }
    

    
    #Flows
    # Adjusted to account for VEp, which reduces the impact of AVEi since it, in essence, reduces fractionSymptomatic
    forceOfInfection <- beta / populationFractions / (fractionSymptomatic + relativeInfectivityAsymptomatic - fractionSymptomatic * relativeInfectivityAsymptomatic) * 
      as.vector((contactMatrix %*% (
          (A + I) * (1 - fractionSymptomatic) * relativeInfectivityAsymptomatic + 
            fractionSymptomatic * (A + (Av*(1 - VEi1) * (1 - VEp1)) + (Avb*(1 - VEi2) * (1 - VEp2))) + 
            (1 - AVEi.eff) * fractionSymptomatic * (I + (Iv * (1 - VEi1) * (1 - VEp1)) + (Ivb * (1 - VEi2) * (1 - VEp2))) + 
            (Av + Iv) * relativeInfectivityAsymptomatic * (1 - VEi1) * (1 - fractionSymptomatic + fractionSymptomatic * VEp1) +
            (Avb + Ivb) * relativeInfectivityAsymptomatic * (1 - VEi2) * (1 - fractionSymptomatic + fractionSymptomatic * VEp2)
          )))



    
    # if(any(S <0)){stop("S < 0")}
    # if(any(Sv <0)){stop(paste("Sv < 0, t =",t))}
    # if(any(Svb <0)){stop("Svb < 0")}

    #if(any(E <0)){stop(paste("E < 0, t = ",t))}
    #if(any(Ev <0)){stop(paste("Ev < 0, t =",t))}
    #if(any(Evb <0)){stop("Evb < 0")}

    # if(any(A <0)){stop("A < 0")}
    # if(any(Av <0)){stop("Av < 0")}
    # if(any(Avb <0)){stop("Avb < 0")}

    # if(any(I <0)){stop("I < 0")}
    # if(any(Iv <0)){stop("Iv < 0")}
    # if(any(Ivb <0)){stop("Ivb < 0")}

    
    
    if(any(forceOfInfection<0)){stop("forceOfInfection < 0 , t = ", t)}
    forceOfInfection=ifelse(forceOfInfection<0,0,forceOfInfection)
    
    S_to_E <- S * forceOfInfection
    E_to_A <- lambda1 * E
    A_to_I <- lambda2 * A
    I_to_R <- gamma * I

    S_to_Sv <-  ifelse(isVaccinatingPrimeByAge, vaccinationRatePrimeByAge * S / EligibleForDose1, 0)
    Sv_to_Ev <- Sv * (1 - VEs1) * forceOfInfection
    Ev_to_Av <- lambda1 * Ev
    Av_to_Iv <- lambda2 * Av
    Iv_to_Rv <- gamma * Iv

    Sv_to_Svb = ifelse(isVaccinatingBoostByAge, vaccinationRateBoostByAge *  Sv / EligibleForDose2, 0)
    Svb_to_Evb <- Svb * (1 - VEs2) * forceOfInfection
    Evb_to_Avb <- lambda1 * Evb
    Avb_to_Ivb = lambda2 * Avb
    Ivb_to_Rvb <- gamma * Ivb
    
    A_and_Av_to_Symp = lambda2 * (A+Av+Avb) * fractionSymptomatic *caseHospitalizationRatio
    Symp_to_Hosp = (1/DelayOnsetToHosp) * Symp 

    Hosp_Discharge = Hosp/hospDuration

    if(any(S_to_Sv<0)){stop("any(S_to_Sv<0)")}
    if(any(Sv_to_Ev<0)){stop("any(Sv_to_Ev<0)")}
    if(any(Ev_to_Av<0)){stop("any(Ev_to_Av<0)")}
    if(any(Av_to_Iv<0)){stop("any(Av_to_Iv<0)")}
    if(any(Iv_to_Rv<0)){stop("any(Iv_to_Rv<0)")}

    if(any(Sv_to_Svb<0)){stop("any(Sv_to_Svb<0)")}
    if(any(Svb_to_Evb<0)){stop("any(Svb_to_Evb<0)")}
    if(any(Evb_to_Avb<0)){stop("any(Evb_to_Avb<0)")}
    if(any(Avb_to_Ivb<0)){stop("any(Avb_to_Ivb<0)")}
    if(any(Ivb_to_Rvb<0)){stop("any(Ivb_to_Rvb<0)")}
    

    
    #Derivatives
    #Non-vaccinated compartments
    dS <- -S_to_E - S_to_Sv
    dE <- S_to_E  - E_to_A
    dA <- E_to_A - A_to_I
    dI <- A_to_I - I_to_R
    dR <- I_to_R
    #Vaccinated compartments
    dSv <- S_to_Sv -Sv_to_Ev - Sv_to_Svb
    dEv <- Sv_to_Ev - Ev_to_Av   
    dAv <- Ev_to_Av - Av_to_Iv
    dIv <- Av_to_Iv - Iv_to_Rv
    dRv <- Iv_to_Rv
    
    dSvb <- Sv_to_Svb -  Svb_to_Evb
    dEvb <- Svb_to_Evb - Evb_to_Avb 
    dAvb <- Evb_to_Avb - Avb_to_Ivb
    dIvb <- Avb_to_Ivb - Ivb_to_Rvb
    dRvb <- Ivb_to_Rvb
    #Auxiliary vaccinated compartment
    dV <- ifelse(isVaccinatingPrimeByAge, vaccinationRatePrimeByAge, 0)
    #Auxiliary vaccinated compartment - people with a boosting dose
    dVb <- ifelse(isVaccinatingBoostByAge, vaccinationRateBoostByAge, 0)
    
    if(any(V==0 & dVb>0)){stop("any(V==0 & dVb>0)")}
    #Hospitalized compartment
    dSymp = A_and_Av_to_Symp - Symp_to_Hosp
    dHosp = Symp_to_Hosp - Hosp_Discharge
    dDischarged = Hosp_Discharge
    zeroVec <- 0 * populationFractions
    
    #Return derivative
    return(list(c(dS, dE, dA, dI, dR, dSv, dEv, dAv, dIv, dRv, dSvb, dEvb, dAvb, dIvb, dRvb, dV, dVb, zeroVec,  zeroVec, dSymp, dHosp, dDischarged)))
  })
}

#This is a utility function that reconstructs the model state as a list so that equations can refer to compartments by name
reconstructState.SEAIRV <- function(state) {
  numberOfClasses <- length(state) / 22 #Each of the 11 classes are vectors of the same length
  S  <- state[                       1 :     numberOfClasses ]
  E  <- state[    (numberOfClasses + 1):(2 * numberOfClasses)]
  A  <- state[(2 * numberOfClasses + 1):(3 * numberOfClasses)]
  I  <- state[(3 * numberOfClasses + 1):(4 * numberOfClasses)]
  R  <- state[(4 * numberOfClasses + 1):(5 * numberOfClasses)]
  Sv <- state[(5 * numberOfClasses + 1):(6 * numberOfClasses)]
  Ev <- state[(6 * numberOfClasses + 1):(7 * numberOfClasses)]
  Av <- state[(7 * numberOfClasses + 1):(8 * numberOfClasses)]
  Iv <- state[(8 * numberOfClasses + 1):(9 * numberOfClasses)]
  Rv <- state[(9 * numberOfClasses + 1):(10 * numberOfClasses)]
  Svb <- state[(10 * numberOfClasses + 1):(11 * numberOfClasses)]
  Evb <- state[(11 * numberOfClasses + 1):(12 * numberOfClasses)]
  Avb <- state[(12 * numberOfClasses + 1):(13 * numberOfClasses)]
  Ivb <- state[(13 * numberOfClasses + 1):(14 * numberOfClasses)]
  Rvb <- state[(14 * numberOfClasses + 1):(15 * numberOfClasses)]
  V  <- state[(15 * numberOfClasses + 1):(16 * numberOfClasses)]
  Vb <- state[(16* numberOfClasses + 1):(17 * numberOfClasses)]
  vaccinatingPrime <- state[(17 * numberOfClasses + 1):(18 * numberOfClasses)]
  vaccinatingBoost <- state[(18 * numberOfClasses + 1):(19 * numberOfClasses)]
 Symp = state[(19 * numberOfClasses + 1):(20 * numberOfClasses)]
  Hosp <- state[(20 * numberOfClasses + 1):(21 * numberOfClasses)]
  Discharged = state[(21 * numberOfClasses + 1):(22 * numberOfClasses)]
  return(as.list(environment()))
}

#This function implements seeding infections in the SEAIR+V model
#parameters should define seedInfections, lambda1, lambda2, and gamma
doSeed.SEAIRV <- function(state, parameters) {
  stateList <- reconstructState.SEAIRV(state)
  with(append(stateList, parameters), {
    seedInfectionsFractions <- seedInfections / population
    S <- S - seedInfectionsFractions
    # E, A, I initialized such that A and I have derivatives of 0 at seed
    E <- E + seedInfectionsFractions * (gamma * lambda2) / 
      (lambda1 * lambda2 + gamma * (lambda1 + lambda2))
    A <- A + 0 * (gamma * lambda1) /
      (lambda1 * lambda2 + gamma * (lambda1 + lambda2))
    I <- I + 0 * (lambda1 * lambda2) /
      (lambda1 * lambda2 + gamma * (lambda1 + lambda2))
    
   # Symp <- (fractionSymptomatic*(I + seedInfectionsFractions) * (lambda1 * lambda2)) /
   #   (lambda1 * lambda2 + gamma * (lambda1 + lambda2))
    
    #Return derivative
    return(c(S, E, A, I, R, Sv, Ev, Av, Iv, Rv,  Svb, Evb, Avb, Ivb, Rvb, V, Vb, vaccinatingPrime, vaccinatingBoost, Symp, Hosp, Discharged))
  })
}

#' #' @title Check antiviral treatment inputs
#' #' @description Checks antiviral treatment inputs and computes the antiviral treatment parameters
#' #' @return List of antiviral treatment parameters
#' #' @keywords internal
#' checkInputs.Antiviral <- function(populationFractions, fractionSymptomatic = 0.5,  fractionSeekCare = 0.6,
#'                                   fractionDiagnosedAndPrescribedOutpatient = 0.7, fractionAdhere = 0.8,
#'                                   fractionAdmitted = 1, fractionDiagnosedAndPrescribedInpatient = 1, AVEi = 0, AVEp = 0, ...) {
#'   #Validate parameters
#'   checkBetween0and1(fractionSymptomatic)
#'   checkDimensionsMatch(fractionSymptomatic, populationFractions)
#'   checkBetween0and1(fractionSeekCare)
#'   checkDimensionsMatch(fractionSeekCare, populationFractions)
#'   checkBetween0and1(fractionDiagnosedAndPrescribedOutpatient)
#'   checkDimensionsMatch(fractionDiagnosedAndPrescribedOutpatient, populationFractions)
#'   checkBetween0and1(fractionAdhere)
#'   checkDimensionsMatch(fractionAdhere, populationFractions)
#'   checkBetween0and1(fractionAdmitted)
#'   checkDimensionsMatch(fractionAdmitted, populationFractions)
#'   checkBetween0and1(fractionDiagnosedAndPrescribedInpatient)
#'   checkDimensionsMatch(fractionDiagnosedAndPrescribedInpatient, populationFractions)
#'   checkBetween0and1(AVEi)
#'   checkDimensionsMatch(AVEi, populationFractions)
#'   checkBetween0and1(AVEp)
#'   checkDimensionsMatch(AVEp, populationFractions)
#'   
#'   AVEi.eff <- AVEi * fractionAdhere * fractionDiagnosedAndPrescribedOutpatient * fractionSeekCare * fractionSymptomatic
#'   return(list(fractionSymptomatic = fractionSymptomatic, fractionSeekCare = fractionSeekCare, 
#'               fractionDiagnosedAndPrescribedOutpatient = fractionDiagnosedAndPrescribedOutpatient, fractionAdhere = fractionAdhere,
#'               fractionAdmitted = fractionAdmitted, fractionDiagnosedAndPrescribedInpatient = fractionDiagnosedAndPrescribedInpatient,
#'               AVEi.eff = AVEi.eff, AVEp = AVEp))
#' }
#' 
#' #' @title Check vaccine inputs
#' #' @description Checks vaccine inputs and computes the vaccine parameters
#' #' @return List of vaccine parameters
#' #' @keywords internal
#' checkInputs.Vaccine <- function(population, populationFractions, seedStartDay, simulationLength,
#'                                 vaccineAdministrationRatePerDay = 0, vaccineAvailabilityByDay = 0,
#'                                 vaccineUptakeMultiplier = 1, VEs = 0, VEi = 0, VEp = 0, 
#'                                 vaccineEfficacyDelay = 7, ...) {
#'   #Validate vaccine parameters
#'   #vaccineAdministrationRatePerDay
#'   checkNonNegativeNumber(vaccineAdministrationRatePerDay)
#'   #vaccineAvailabilityByDay
#'   checkNonNegative(vaccineAvailabilityByDay)
#'   #vaccineUptakeMultiplier
#'   checkNonNegative(vaccineUptakeMultiplier)
#'   checkDimensionsMatch(vaccineUptakeMultiplier, populationFractions)
#'   #VEs
#'   checkBetween0and1(VEs)
#'   checkDimensionsMatch(VEs, populationFractions)
#'   #VEi
#'   checkBetween0and1(VEi)
#'   checkDimensionsMatch(VEi, populationFractions)
#'   #VEp
#'   checkBetween0and1(VEp)
#'   checkDimensionsMatch(VEp, populationFractions)
#'   #vaccineEfficacyDelay
#'   checkNonNegativeNumber(vaccineEfficacyDelay)
#'   
#'   
#'   #Compute the daily vaccination rate
#'   totalSimulationLength <- seedStartDay + simulationLength
#'   vaccinationRateByDay <- rep(0, totalSimulationLength)
#'   currentVaccineAvailability <- 0
#'   for (i in 1:totalSimulationLength) {
#'     if (i <= length(vaccineAvailabilityByDay)){
#'       currentVaccineAvailability <- currentVaccineAvailability + vaccineAvailabilityByDay[i]
#'     }
#'     vaccinationRateByDay[i] <- min(vaccineAdministrationRatePerDay, currentVaccineAvailability)
#'     currentVaccineAvailability <- currentVaccineAvailability - vaccinationRateByDay[i]
#'   }
#'   vaccinationRateByDay <- vaccinationRateByDay / population #Normalize
#'   
#'   #Define vaccination rate function
#'   vaccinationRate <- function(t) {
#'     if ((t < vaccineEfficacyDelay) || (t >= totalSimulationLength + vaccineEfficacyDelay)) {
#'       return(0)
#'     } else {
#'       return(vaccinationRateByDay[floor(t - vaccineEfficacyDelay + 1)])
#'     }
#'   }
#'   
#'   #Compute the vaccination rate age multiplier
#'   vaccinationRateAgeMultiplier <- vaccineUptakeMultiplier * populationFractions
#'   totalMultiplier <- sum(vaccinationRateAgeMultiplier)
#'   if (totalMultiplier > 0) {
#'     vaccinationRateAgeMultiplier <- vaccinationRateAgeMultiplier / totalMultiplier
#'   } else {
#'     warning("vaccineUptakeMultiplier prevents vaccination from occurring.", call. = FALSE)
#'   }
#'   
#'   #Return the parameters
#'   return(list(vaccinationRate = vaccinationRate, vaccinationRateAgeMultiplier = vaccinationRateAgeMultiplier,
#'               VEs = VEs, VEi = VEi, VEp = VEp, vaccineEfficacyDelay = vaccineEfficacyDelay))
#' }