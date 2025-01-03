# A function containing calls to the functions necessary to initialize the model parameters

parameterInitialisation <- function(par, vegData_stage1, vegData_stage2){
  
  # The bulk density profile is calculated
  par@bulkDensityProfile <- calculateBulkDensity(par)
  
  # Depth profiles of Km values
  par <- kmDepthProf(par)
  
  # Depth profiles of CO2 diffusivity
  par <- calculate_CO2diff(par)
  
  # The depth profiles of root carbon inputs
  out <- rootInputs(par, vegData_stage1, vegData_stage2)
  vegData_stage1 <- out$vegData_stage1
  vegData_stage2 <- out$vegData_stage2
  rm(out)
  
  # The annual root CO2 production depth profile
  out <- annualrootCO2Production(par, vegData_stage1, vegData_stage2)
  vegData_stage1 <- out$vegData_stage1
  vegData_stage2 <- out$vegData_stage2
  rm(out)
  
  # The biodiffusion coefficient depth profile
  par <- bioDiff(par)
  
  # The fraction of soil influenced by the rhizosphere
  out <- rhizosphereVolume(par, vegData_stage1, vegData_stage2)
  vegData_stage1 <- out$vegData_stage1
  vegData_stage2 <- out$vegData_stage2
  rm(out)
  
  # Adsorption and desorption parameters are calculated
  par <- adsParam(par, vegData_stage1, vegData_stage2)
  
  # The isotopic data are loaded
  if(par@calibMode == 0){
    
    par <- loadD13CData(par, vegData_stage1, vegData_stage2)    # Loading the d13C data
    par <- loadD14CData(par)                                    # Loading the d14C data
    out <- constructIsotopeInputs(par, vegData_stage1, vegData_stage2)  # Construct the annual inputs of 12C, 13C and 14C
    isotopicInputs <- out$isotopicInputs
    par <- out$par
    isotopeRatios <- out$isotopeRatios
    rm(out)
    
  } else if(par@calibMode == 1 & par@counter_loadIsotopes == 0){
    
    par <- loadD13CData(par, vegData_stage1, vegData_stage2)    # Loading the d13C data
    par <- loadD14CData(par)                                    # Loading the d14C data
    out <- constructIsotopeInputs(par, vegData_stage1, vegData_stage2)  # Construct the annual inputs of 12C, 13C and 14C
    isotopicInputs <- out$isotopicInputs
    par <- out$par
    isotopeRatios <- out$isotopeRatios
    rm(out)
    
  }
  
  # The depth profiles of 12CO2, 13CO2 and 14CO2 are initialized
  par <- initialCO2profile(par, vegData_stage1)
  
  # The state variables are initialized
  x_init <- initialization(par, vegData_stage1, isotopeRatios, isotopicInputs)
  
  if(par@calibMode == 0){
    
    output <- list("par" = par, 
                   "vegData_stage1" = vegData_stage1, 
                   "vegData_stage2" = vegData_stage2, 
                   "isotopicInputs" = isotopicInputs, 
                   "isotopeRatios" = isotopeRatios, 
                   "x_init" = x_init)
    
  } else if(par@calibMode == 1 & par@counter_loadIsotopes == 0){
    
    output <- list("par" = par, 
                   "vegData_stage1" = vegData_stage1, 
                   "vegData_stage2" = vegData_stage2, 
                   "isotopicInputs" = isotopicInputs, 
                   "isotopeRatios" = isotopeRatios, 
                   "x_init" = x_init)
    
  } else if(par@calibMode == 1 & par@counter_loadIsotopes == 1){
    
    output <- list("par" = par)
    
  }

  return(output)
  
}