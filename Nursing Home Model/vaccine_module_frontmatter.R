## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------
# knitr::opts_chunk$set(echo = FALSE)



## ----load required packages and define helper functions, warning=FALSE, message=FALSE--------------------------------------------------------------------------------------
# try(setwd(wd<-"C:\\Users\\prabasaj\\OneDrive - CDC\\corona"))
# try(setwd(wd<-"C:\\Users\\vig5\\OneDrive - CDC\\corona"))
options(stringsAsFactors=FALSE)
library(rio)
library(deSolve)

# Country specific social contact rate data (suppl to Prem et al)
polymod = import("contact_matrices_152_countries/MUestimates_all_locations_2.xlsx",
  which="United States of America", col_names=FALSE)
polymod.home = import("contact_matrices_152_countries/MUestimates_home_2.xlsx",
  which="United States of America", col_names=FALSE)
polymod.school = import("contact_matrices_152_countries/MUestimates_school_2.xlsx",
  which="United States of America", col_names=FALSE)
polymod.work = import("contact_matrices_152_countries/MUestimates_work_2.xlsx",
  which="United States of America", col_names=FALSE)
polymod.other = import("contact_matrices_152_countries/MUestimates_other_locations_2.xlsx",
  which="United States of America", col_names=FALSE)

## Load county data for capacity analysis
# Healthcare capacity and use data:
load("county_capacity_hrsa_aha.RData")
hosp_capacity = read.delim("hosp_capacity.txt")
hrr_hosp = import("BEDS_ICU_POPULATION_HRR_2017_2018.xlsx")
state_hosp = import("BEDS_ICU_VENTS_STATE_formatted.xlsx")
names(hrr_hosp) = toupper(names(hrr_hosp))

# ACS 2018 5yr population by age (and sex)
acspop = read.delim("ACStable_B01001_50_2018_5.txt")
acspop$GEOID = substr(10^5 + acspop$STATE*1000 + acspop$COUNTY, 2, 6)
ages.in.col = lapply(strsplit(gsub("(^[[:alpha:].]+)|([.]years.*$)|(ZCTA5)", "", names(acspop)), "[^[:digit:]]+"), function(ll) {
  res = ifelse(length(ll)==0 | any(is.na(as.numeric(ll))), NA, list(as.numeric(ll[1]):as.numeric(tail(ll,1))))[[1]]
  if (identical(res, as.integer(5))) res = 0:4
  if (identical(res[[1]], as.integer(85))) res = 85:100
  res
})

get_pops16_ACS = function(geoid = "13089") { # DeFault DeKalb GA
  acs2poly = lapply(2:length(polymod.cuts), function(ii) polymod.cuts[ii-1]:min(polymod.cuts[ii]-1, 100))
  age_cols = lapply(acs2poly, function(ll16) sapply(ages.in.col, function(ll) length(intersect(ll16, ll))>0))
  if (any(rowSums(data.frame(age_cols))>1)) print("Warning: supplied cuts do not respect ACS age cuts in get_pops16_ACS")
  myset = subset(acspop, GEOID %in% geoid | identical(geoid, "US")) 
  sapply(age_cols, function(cc) sum(myset[,cc]))
}

get_popsUNI_ACS = function(geoid = "13089") { # DeFault DeKalb GA
  age_cols = list(HS=sapply(ages.in.col, function(ll) any(0:17 %in% ll)), 
                  UNI=sapply(ages.in.col, function(ll) any(18:21 %in% ll)),
                  OLD=sapply(ages.in.col, function(ll) any(22:100 %in% ll)))
  age_strata = names(age_cols)
  myset = subset(acspop, GEOID %in% geoid | identical(geoid, "US")) 
  pops = sapply(age_cols, function(cc) sum(myset[,cc]))
  names(pops) = age_strata
  pops
}

trans_prob_to_doubling_time = function(trans_per_cont, Pars) {
  # REM: beta_ij S_j I_i in naive population -> t(trans_per_cont*poly.US) %*% I
  nr = length(Pars$age_frac)
  b11 = t(trans_per_cont*Pars$inf_E_over_I*poly.US) - diag(1/Pars$tau_E, nr)
  b12 = t(trans_per_cont*poly.US)
  b21 = diag(1/Pars$tau_E, nr)
  b22 = -diag(1/Pars$tau_I, nr)
  b = rbind(cbind(b11, b12), cbind(b21, b22))
  res = log(2)/max(eigen(b)$values)
  ifelse(res<0, Inf, res)
}

# Doubling time and principal eigenvector from next generation matrix
next_gen_analysis = function(Pars, verbose=FALSE) {
  with(Pars, {
    nr = nrow(beta_I)
    b11 = age_frac*t(beta_E) - diag(1/tau_E, nr)
    b12 = age_frac*t(beta_I)
    b21 = diag(1/tau_E, nr)
    b22 = -diag(1/tau_I, nr)
    b = rbind(cbind(b11, b12), cbind(b21, b22))
    ev = eigen(b)
    mval = max(Re(ev$values))
    mvec = ev$vectors[, which.max(Re(ev$values))]
    check.sign = (sum(abs(sign(Re(mvec)))) == abs(sum(sign(Re(mvec))))) # Check if max vector is all of same sign or zero
    if (!check.sign & verbose) print("Anomalous eigenvector of largest eigenvalue detected")
    list(doubling_time=log(2)/max(0, Re(eigen(b)$values)), eval=ev$values[which.max(Re(ev$values))], evec=mvec)
  })
}

pop_adjust_pars = function(Pars, geoid="US", norm_100K=FALSE, betas=FALSE) {
  # Boiler plate updates
  Pars$max_day = 365 # Simulation period
  Pars$tau_I = with(Pars,  tau_E * inf_E_over_I * (1-prop_E_trans)/prop_E_trans)
  Pars$rate_hosp = with(Pars, (symp_case_hosp_frac/(1 - symp_case_hosp_frac))/tau_I) # get_hosp applies rate to frac_sym, 1/(1-.) correction 2020-03-31
  if (!is.null(Pars$generation_time) & is.null(Pars$doubling_time)) {
    Pars$doubling_time = with(Pars, log(2) * generation_time/log(R0))
  }

  Pars$GEOID = geoid
  Pars$population = get_pops16_ACS(Pars$GEOID) # Alt: get_pops_ACS(Pars$GEOID)
  Pars$age_frac = prop.table(Pars$population) # Proportion in age group

  poly.US <<- polymod.US
  Pars$prob_trans_per_contact = uniroot(function(x) trans_prob_to_doubling_time(x, Pars)-Pars$doubling_time, c(0, 1))$root


  if (norm_100K) Pars$population = 10^5 * Pars$age_frac # Normalize base population to 100K
  
  if (betas) {
    # Default (baseline) matrices, modified for different scenarios:
    Pars$polymod = with(Pars, t(t(polyscale.US) * age_frac)) # Density-rescaling
    # Remember that beta_ij is for a 100% population of naive j
    Pars$beta_I = with(Pars, prob_trans_per_contact * polyscale.US)
    Pars$beta_E = with(Pars, inf_E_over_I * beta_I)
  }

  Pars
}

get_hosp_capacity = function(geoid) {
  hrsa_beds = sum(subset(hrsa.sm, f00002 %in% geoid | identical(geoid, "US"))$f0892117, na.rm=TRUE)
  hrsa_beds_in_use = sum(subset(hrsa.sm, f00002 %in% geoid | identical(geoid, "US"))$f0954517, na.rm=TRUE)/365
  hrsa_icu_beds = with(subset(hrsa.sm, f00002 %in% geoid | identical(geoid, "US")), 
    sum(f0913917+f0913317+f0914517+f0916317+f0912717+f0910317+f1330917, na.rm=TRUE))
  aha_beds = with(subset(aha, GEOID %in% geoid | identical(geoid, "US")), sum(Beds, na.rm=TRUE))
  aha_beds_in_use = with(subset(aha, GEOID %in% geoid | identical(geoid, "US")), sum(Beds*Occupancy, na.rm=TRUE))
  ltach_beds =  with(subset(hrsa.sm, f00002 %in% geoid | identical(geoid, "US")), 
    sum(ltach_beds, na.rm=TRUE))
  vsnf_beds =  with(subset(hrsa.sm, f00002 %in% geoid | identical(geoid, "US")), 
    sum(vsnf_beds, na.rm=TRUE))
  list(hrsa_beds=hrsa_beds, hrsa_beds_in_use=hrsa_beds_in_use, hrsa_icu_beds=hrsa_icu_beds,
       aha_beds=aha_beds, aha_beds_in_use=aha_beds_in_use,
       vsnf_beds=vsnf_beds, ltach_beds=ltach_beds)
}

# Post processing: hospitalized
get_hosp = function(out_state, Pars) {
    symp = cbind(Day=1:dim(out_state)[1], data.frame(t(apply(out_state[,, "I"], 1, function(rr) Pars$frac_symp * rr))))
    hosp = symp[1,]
    hosp[, -1] = (1 - exp(-Pars$rate_hosp*hosp$Day[1])) * hosp[, -1] # Zero on day 0
    hosp.symp = hosp # Hospitalized symptomatics; added 2020-03-31
    daily_adm = c()
    for (rr in 2:nrow(symp)) {
      adm = (symp[rr-1, -1] - hosp.symp[rr-1, -1]) * (1 - exp(-Pars$rate_hosp*diff(symp$Day[rr-1:0])))
      not.dschg = hosp[rr-1, -1] * exp(-Pars$rate_dschg*diff(symp$Day[rr-1:0]))
      not.dschg.symp = not.dschg * exp(-diff(symp$Day[rr-1:0])/Pars$tau_I) # Not discharged and still symptomatic
      hosp.symp = rbind(hosp.symp, c(Day=symp$Day[rr], not.dschg.symp + adm))
      hosp = rbind(hosp, c(Day=symp$Day[rr], not.dschg + adm))
      daily_adm = rbind(daily_adm, adm)
    }
    list(hosp, daily_adm)
}



## ----defining the main ODE solving wrappers, warning=FALSE, message=FALSE, echo=TRUE---------------------------------------------------------------------------------------

age_SEIR = function(Time, State, Pars) {
  # State is 4N-element vector (N-age-group X 4-SEIR)
  # Each element is number of persons/total population (sum(State)=1)
  # ODEs set up to work on N X 4 matrix
  nr = nrow(Pars$beta_I) # Gettin the number of age groups, N
  State = matrix(State, nrow=nr, dimnames=list(paste0("AGE", 1:nr), c("S", "E", "I", "R"))) 
  dState = State * 0
  beta_mult = 1
  if (!is.null(Pars$seasonality)) beta_mult = with(Pars$seasonality, trans_prob_multiplier(Time))
  dState[, "S"] = beta_mult * with(Pars, - (State[, "E"] %*% beta_E) * State[, "S"] - (State[, "I"] %*% beta_I) * State[, "S"])
  dState[, "E"] = with(Pars, - dState[, "S"] - State[, "E"]/tau_E)
  dState[, "I"] = with(Pars, State[, "E"]/tau_E - State[, "I"]/tau_I)
  dState[, "R"] = with(Pars, State[, "I"]/tau_I)
  list(dState)
}

### Converting from scenarios to betas

runODE = function(Pars) {
  # Naive state: all susceptible
  state_naive = cbind(S=Pars$age_frac, E=0*Pars$age_frac, I=0*Pars$age_frac, R=0*Pars$age_frac)
  state_naive = state_naive/sum(state_naive) # For good measure

  # Initial state may be supplied through Pars
  if (is.null(Pars$init_state)) {
    # Initial state: seeded by default by a 35-39 year old (why not?)
    state = state_naive
    midgp = round(length(Pars$age_frac)/2)
    state[midgp, c("S", "E")] = state[midgp, "S"] * (c(1,0) + c(-1, 1)/Pars$population[midgp]) # 1 of 35-39 age group exposed
  } else {
    state = Pars$init_state
  }

  Pars = scenario2beta(Pars)$Pars

  if (is.null(Pars$init_day)) Pars$init_day = 1
  if (is.null(Pars$step_day)) Pars$step_day = 1 # step size, for output

  # Run ODE solver
  days = seq(Pars$init_day, Pars$max_day, by=Pars$step_day)
  event_days = Pars$vax_pars$vax_day + Pars$vax_pars$t_d + c(0, Pars$vax_pars$t_12) # Event = immunogenicity achieved
  event_days = intersect(event_days, days)
  if (length(event_days) > 0) {
    out = ode(age_SEIR, y = c(state), times = days, parms=Pars, events=list(func=Pars$vaccinate, time=event_days))
  } else {
    out = ode(age_SEIR, y = c(state), times = days, parms=Pars)
  }
  
  # Format output to age-stratified 3D array
  out_state = array(out[,-1], dim=c(nrow(out), dim(state_naive)))
  dimnames(out_state) = list(time=paste("Day", out[,1]), age=dimnames(state_naive)[[1]], state=dimnames(state_naive)[[2]])

  # Post-processing: incident cases (daily count) of I and E
  inc_I = apply(out_state, 1, function(s) s[,"E"]/Pars$tau_E) * sum(Pars$population)
  inc_E = apply(out_state, 1, function(s) {
    (s[, "E"] %*% Pars$beta_E) * s[, "S"] + (s[, "I"] %*% Pars$beta_I) * s[, "S"]
    }) * sum(Pars$population)
  peak_incidence = max(colSums(matrix(inc_I, nrow=dim(out_state)[2]))) # Ensuring matrix for single age-group models
  peak_incidence_date = days[which.max(colSums(matrix(inc_I, nrow=dim(out_state)[2])))]

  list(out_state=out_state, inc_E=inc_E, inc_I=inc_I, peak_incidence=peak_incidence, peak_incidence_date=peak_incidence_date, 
    Pars=Pars)
}



## ----Symmetrize and population adjust polymod, warning=FALSE, message=FALSE------------------------------------------------------------------------------------------------
#######################################################
polymod.cuts = c(0:15 * 5, Inf)
US16 = get_pops16_ACS("US") # Use ACS data
# Regularize the darned polymod matrix for the national population distribution
# Arithmetic mean, to ensure additivity and consistency with flu-models approach
polymod.US = ((US16 * polymod) + t(US16 * polymod))/(2 * US16)
polymod.school.US = ((US16 * polymod.school) + t(US16 * polymod.school))/(2 * US16)
polymod.work.US = ((US16 * polymod.work) + t(US16 * polymod.work))/(2 * US16)
polymod.home.US = ((US16 * polymod.home) + t(US16 * polymod.home))/(2 * US16)
polymod.other.US = ((US16 * polymod.other) + t(US16 * polymod.other))/(2 * US16)

# Separating out PK12 from university in school polymod (sweeping - but reasonable? - assumptions)
# Estimating polymod from polymod.school and population fractions in university
# Assume that p15to17 * school interactions 15-19 are at HS, remainder at university
# Assume all school interactions 20-24 are at universities
# Assume that only school interactions with 18-22 are fractionated out for the remainder of age groups as at university
USUNI = get_popsUNI_ACS("US")
p15to17 = (USUNI[1] - sum(US16[1:3]))/US16[4]
p22to24 = (sum(US16[1:5]) - sum(USUNI[1:2]))/US16[5]
polymod.UNI = c(0, 0, 0, 1-p15to17, p22to24, rep(0,11)) * US16 * as.matrix(polymod.school.US)
# All 20-24 school interactions are at university (assumption):
polymod.UNI[5,5] = US16[5] * polymod.school.US[5,5]
# Symmetrize number of 15-19 to 20-24 intercations:
polymod.UNI[4,5] = (polymod.UNI[5,4] + polymod.UNI[4,5])/2
polymod.UNI[5,4] = polymod.UNI[4,5]
# Symmetrize the rest:
polymod.UNI[,4:5] = t(polymod.UNI[4:5,])
# Normalize to population
polymod.UNI = polymod.UNI/US16
polymod.PK12 = polymod.school.US - polymod.UNI

## Scaled polymods -- very important!
polyscale.US = t(t(polymod.US)/prop.table(US16)) # m_ij/P_j
polyscale.home.US = t(t(polymod.home.US)/prop.table(US16)) # m_ij/P_j
polyscale.work.US = t(t(polymod.work.US)/prop.table(US16)) # m_ij/P_j
polyscale.school.US = t(t(polymod.school.US)/prop.table(US16)) # m_ij/P_j
polyscale.other.US = t(t(polymod.other.US)/prop.table(US16)) # m_ij/P_j
polyscale.UNI = t(t(polymod.UNI)/prop.table(US16)) # m_ij/P_j
polyscale.PK12 = t(t(polymod.PK12)/prop.table(US16)) # m_ij/P_j



## ----set up parameters, warning=FALSE, message=FALSE, echo=TRUE------------------------------------------------------------------------------------------------------------
#######################################################
# Set up parameters
# Updated 2020-04-02
pars_scenarios = readRDS("Proposed Outbreak Scenarios COVID 19 current.rds")

# Set population adjusted parameters
pars = pars_scenarios$Best
pars = pop_adjust_pars(pars, "US", norm_100K=TRUE)



## ----Load standard scenarios to compare, warning=FALSE, message=FALSE------------------------------------------------------------------------------------------------------
##########################################################################################
# Some "standard" social distancing scenarios and a helper function

# Each scenario is implemented 
scenario2beta = function(Pars, mult=FALSE) {
  ## Adjusts for both time reallocation among contexts and distancing in each context
  ##  by adjusting age-specific weights on context-specific polymods while maintaining reciprocity
  age_map = list(child=1:4, adult=5:13, `65+`=14:16, all=1:16)
  # Estimate proportion of interaction-time spent in each setting (column) by age group (row); each row sums to 1 
  polymod3D = array(unlist(c(polymod.home.US, polymod.work.US, polymod.school.US, polymod.other.US)), dim=c(dim(polymod.home.US), 4))
  dimnames(polymod3D)[[3]] = c("home","work","school", "other")
  Q_ia = zapsmall(prop.table(apply(polymod3D, 3, rowSums),1)) # Proportion of all interactions
  # Alt: zapsmall(prop.table(apply(polymod3D, 3, diag),1)) # Proportion of interactions within age-group; noisier?

  # Build the two polymod multiplier matrices Pt and qoverQ
  Pt = polymod3D * 0 + 1
  qoverQ = Q_ia * 0 + 1
  sc = Pars$scenario
  for (rr in 1:nrow(sc)) {
    if (sc$`Coefficient-type`[rr] == "qoverQ") {
      qoverQ[age_map[[sc$`Age-group`[rr]]], sc$Context[rr]] = sc$Value[rr]
    }
    if (sc$`Coefficient-type`[rr] == "Pt") {
      Pt[age_map[[sc$`Age-group`[rr]]], age_map[[sc$`Age-group2`[rr]]], sc$Context[rr]] = sc$Value[rr]
      Pt[age_map[[sc$`Age-group2`[rr]]], age_map[[sc$`Age-group`[rr]]], sc$Context[rr]] = sc$Value[rr]
    }
  }
  q_ia = prop.table(qoverQ * Q_ia, 1) # Normalize
  qoverQ = q_ia/Q_ia # Recalculate ratio; will not, in general stay the same
  qoverQ = sqrt(array(apply(qoverQ, 2, function(vv) outer(vv, vv)), dim=dim(polymod3D))) # Expand out to sqrt(q_iaq_ja/Q_ia_Q_ja)
  dimnames(qoverQ)[[3]] = dimnames(polymod3D)[[3]]
  # Convert to actual probabilities
  Pt = Pars$prob_trans_per_contact * Pt
  
  # Build betas from definition
  Pars$beta_E = Pars$inf_E_over_I * (
    Pt[,,"home"] * qoverQ[,, "home"] * polyscale.home.US +
    Pt[,,"school"] * qoverQ[,, "school"] * polyscale.school.US +
    Pt[,,"work"] * qoverQ[,, "work"] * polyscale.work.US +
    Pt[,,"other"] * qoverQ[,, "other"] * polyscale.other.US)       
  Pars$beta_I = (
    Pt[,,"home"] * qoverQ[,, "home"] * polyscale.home.US +
    Pt[,,"school"] * qoverQ[,, "school"] * polyscale.school.US +
    Pt[,,"work"] * qoverQ[,, "work"] * polyscale.work.US +
    Pt[,,"other"] * qoverQ[,, "other"] * polyscale.other.US)    

  list(Pars=Pars, Pt=Pt, qoverQ=qoverQ, q_ia=q_ia)
}

social_distancing_scenarios = import("CM_NPI_coefficients.xlsx")
reallocate_time = FALSE
if (!reallocate_time) social_distancing_scenarios$`Coefficient-type` = "Pt"



## ----social distancing table-----------------------------------------------------------------------------------------------------------------------------------------------
# knitr::kable(social_distancing_scenarios)



## ----large bundle of scenarious by trigger and duration, warning=FALSE, message=FALSE--------------------------------------------------------------------------------------
###########################################################################
# Implement seasonality
pars$seasonality = list(t_winter=0, # Origin for sinusoid peak (0 = start of outbreak)
                        amplitude=0.3, # (peak-valley)/peak; (2-1.4)/2 from Lipsitch et al (speculative)
                        time_period=365, # year, for seasonality, but modifiable for other effects
                        trans_prob_multiplier = function(time) {
                          amp = get("amplitude", parent.frame())
                          t_w = get("t_winter", parent.frame())
                          t_p = get("time_period", parent.frame())
                          (1 - abs(amp)/2) + (abs(amp)/2) * cos(2*pi * (time-t_w) / t_p)}
)
   
##########################################################################
# Set vaccine parameters and function
pars$vax_pars = list(eff1=0.2, # First dose effectiveness
                     eff2=0.95, # Second dose effectiveness
                     t_12=28,  # Time (days) between doses
                     t_d=14,   # Time (days) from dose to corresponding immunogenicity
                     sero_test=TRUE, # If serological testing before vaccination
                     vax_day=228, # Day of first dose (single day, for now)
                     vax_doses=rep(7500/16E5, 16) # Doses, by age-group, normalized to population
)

pars$vaccinate = function(t, y, parms) { # Vaccination events$func
  y = matrix(y, ncol=4)
  dimnames(y)[[2]] = c("S","E","I","R")
  vac_frac =  
    parms$vax_pars$vax_doses * y[,"S"] * ifelse(t %in% (parms$vax_pars$vax_day + parms$vax_pars$t_d), parms$vax_pars$eff1, parms$vax_pars$eff2 - parms$vax_pars$eff1) /
    (y[,"S"] + y[,"E"] + (!parms$vax_pars$sero_test) * (1-parms$frac_symp) * (y[,"I"] + y[,"R"])) 
  vac_frac = pmin(vac_frac, y[,"S"])
  y[,"S"] = y[,"S"] - vac_frac
  y[,"R"] = y[,"R"] + vac_frac
  as.vector(y)
}

pandemic_plan = htmltab::htmltab("https://www.cdc.gov/flu/pandemic-resources/national-strategy/planning-guidance/appendix-a.html")
pandemic_plan$`Estimated Number in Group` = as.numeric(gsub(",", "", pandemic_plan$`Estimated Number in Group`))
pandemic_pop5 = sum(pandemic_plan$`Estimated Number in Group`[1:5])/sum(pandemic_plan$`Estimated Number in Group`)
pandemic_vaxxed = 25E6/sum(pandemic_plan$`Estimated Number in Group`)

vax_dose_schemes = list()
vax_dose_schemes[[1]] = pars$age_frac * 0 # Unvaxxed
vax_old = TRUE; vax_child = TRUE
vax_dose_schemes[[2]] = prop.table(pars$age_frac * c(rep(0,4), rep(1,9), rep(0,3))) * pandemic_pop5  + 
                  prop.table(pars$age_frac * c(rep(vax_child,4), rep(0,9), rep(vax_old,3))) * (pandemic_vaxxed - pandemic_pop5)
vax_old = FALSE; vax_child = TRUE
vax_dose_schemes[[3]] = prop.table(pars$age_frac * c(rep(0,4), rep(1,9), rep(0,3))) * pandemic_pop5  + 
                  prop.table(pars$age_frac * c(rep(vax_child,4), rep(0,9), rep(vax_old,3))) * (pandemic_vaxxed - pandemic_pop5)
vax_old = TRUE; vax_child = FALSE
vax_dose_schemes[[4]] = prop.table(pars$age_frac * c(rep(0,4), rep(1,9), rep(0,3))) * pandemic_pop5  + 
                  prop.table(pars$age_frac * c(rep(vax_child,4), rep(0,9), rep(vax_old,3))) * (pandemic_vaxxed - pandemic_pop5)


pars$scenario = subset(social_distancing_scenarios, Scenario=="Baseline")


vax_run = function(input) {
  target = input$target
  pars$seasonality$amplitude = input$amp # The "Lipsitch default" is 0.3
  pars$max_day = 730 # 365 usually; 730 To ensure seasonality plays out

  pars$vax_pars$vax_day = Inf
  pars = suppressWarnings(within(pars, rm(init_state, init_day)))
  pars$scenario = subset(social_distancing_scenarios, Scenario=="Baseline")
  baserun = runODE(pars)
  trig_dt = which(rowSums(baserun$out_state[,,"I"] + baserun$out_state[,,"R"])>target)[1]

  # closure_durations = c(6, 8, 12, 16) * 7 
  # vax_days = c(170, 200) + trig_dt
  # vax_schemes = seq(vax_dose_schemes) 
  # all_sc = setdiff(unique(social_distancing_scenarios$Scenario), "Lenient")

  pps = input$pps
  closure_duration = max(input$closure_duration * 7, 1) # Prevents broken indexing for 0 day closure
  vax_scheme = input$vax_scheme
  vax_day = input$vax_day - as.Date("2020-03-15") + trig_dt

  post_pause_scenario = subset(social_distancing_scenarios, Scenario==pps)
  pars$vax_pars$vax_doses = vax_dose_schemes[[as.numeric(vax_scheme)]]
  pars$vax_pars$vax_day = vax_day
  
  if (is.null(pars$pause_scenario)) pars$pause_scenario = "Surge" 
  pars$scenario = subset(social_distancing_scenarios, Scenario==pars$pause_scenario)
  pars$init_state = c(baserun$out_state[trig_dt,,])
  pars$init_day = trig_dt
  out_state2 = runODE(pars)$out_state

  pars$scenario = post_pause_scenario
  pars$init_state = c(out_state2[closure_duration+1,,])
  pars$init_day = trig_dt + closure_duration
  out_state3 = runODE(pars)$out_state
  pars$vax_pars$vax_doses = 0 * pars$vax_pars$vax_doses
  out_state3_novax = runODE(pars)$out_state

  out_state = baserun$out_state
  out_state[trig_dt - 1 + 1:closure_duration,,] = out_state2[1:closure_duration,,]
  out_state[(trig_dt + closure_duration):pars$max_day,,] = out_state3
  out_state_novax = out_state
  out_state_novax[(trig_dt + closure_duration):pars$max_day,,] = out_state3_novax

  list(baserun=baserun$out_state, out_state=out_state, out_state_novax=out_state_novax, trig_dt=trig_dt, input=input)
}

