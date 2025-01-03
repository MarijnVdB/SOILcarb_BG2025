# This script is used to load the 13C input data

loadD13CData <- function(par, vegData_stage1, vegData_stage2){
  
  #If the 'Terrestrial Suess effect is included, this data is loaded
  if(par@Include_Suess_effect == 1){
    
    # The data is loaded
    allData_d13C <- read.csv("SuessEffectD13CDifference.csv") # Historic atmospheric d13CO2
    # allData_d13C <- read.csv("SuessEffectD13CDifference_const_d13CO2.csv") # Historic atmospheric d13CO2
    allData_CO2 <- read.csv("Historic_CO2_data.csv")          # Historic atmospheric [CO2]
    
    # The last year for which Suess differences are available
    lastYear_d13C <- max(allData_d13C$Year)
    lastYear_CO2Conc <- max(allData_CO2$Year)
    
    # If the user provides a last 'real' simlation year that is more recent
    # than the last year for which Suess data is available, an error is
    # generated and the program is stopped 
    if(par@lastRealSimulationYear > lastYear_d13C){
      stop("You provided a lastRealSimulationYear that is more recent then the last year for which Suess data is available!")
    }
    if(par@lastRealSimulationYear > lastYear_CO2Conc){
      stop("You provided a lastRealSimulationYear that is more recent then the last year for which atmospheric CO2 data data is available!")
    }
    
    #The range of atmospheric d13CO2 and [CO2] is adapted for the last 'real' simulation year
    rowNum <- which(allData_d13C$Year == par@lastRealSimulationYear)
    allData_d13C <- allData_d13C[1:rowNum,]
    
    rowNum <- which(allData_CO2$Year == par@lastRealSimulationYear)
    allData_CO2 <- allData_CO2[1:rowNum,]
    
    # Atmospheric d13CO2 throughout the simulation
    d13C_atmosphere <- rep(NA,par@numberOfSimulationYears)
    
    # If more years are simulated than d13CO2 is provided
    if(par@numberOfSimulationYears > nrow(allData_d13C)){
      endOf <- length(d13C_atmosphere)
      d13C_atmosphere[(endOf-nrow(allData_d13C)+1):endOf] <- allData_d13C$d13CO2
      d13C_atmosphere[1:(endOf-nrow(allData_d13C))] <- allData_d13C$d13CO2[1]
    } else{ # If for all simulated years d13CO2 data is provided
      endOf <- nrow(allData_d13C)
      d13C_atmosphere <- allData_d13C$d13CO2[(endOf-par@numberOfSimulationYears+1):endOf]
    }
    
    # Atmospheric [CO2] throughout the simulation
    CO2conc_atmosphere <- rep(NA,par@numberOfSimulationYears)
    # If more years are simulated than [CO2] is provided
    if(par@numberOfSimulationYears > nrow(allData_CO2)){
      endOf <- length(CO2conc_atmosphere)
      CO2conc_atmosphere[(endOf-nrow(allData_CO2)+1):endOf] <- allData_CO2$CO2
      CO2conc_atmosphere[1:(endOf-nrow(allData_CO2))] <- allData_CO2$CO2[1]
    } else{ # If for all simulated years [CO2] data is provided
      endOf <- nrow(allData_CO2)
      CO2conc_atmosphere <- allData_CO2$CO2[(endOf-par@numberOfSimulationYears+1):endOf]
    }
    
    # The D13C is calculated: the difference in vegetation d13C and
    # atmospheric d13CO2; based on Keeling et al. (2017)
    
    d13C_AGveg <- rep(NA,par@numberOfSimulationYears)
    d13C_BGveg <- rep(NA,par@numberOfSimulationYears)
    d13C_rhizodeposits <- rep(NA,par@numberOfSimulationYears)
    D13CSuessDifference_AGveg <- rep(NA,par@numberOfSimulationYears)
    D13CSuessDifference_BGveg <- rep(NA,par@numberOfSimulationYears)
    D13CSuessDifference_rhizodeposits <- rep(NA,par@numberOfSimulationYears)
    
    # The difference in d13C between plant and d13CO2 is corrected for the [CO2],
    # based on Keeling et al. (2017)
    D13CSuessDifference_AGveg[1:par@yearsStage1] <- 
      (tail(d13C_atmosphere,1) - vegData_stage1@AGveg_d13C) - 
      (tail(CO2conc_atmosphere,1) - CO2conc_atmosphere[1:par@yearsStage1])*par@CO2conc_fractionationFactor
    D13CSuessDifference_BGveg[1:par@yearsStage1] <- 
      (tail(d13C_atmosphere,1) - vegData_stage1@BGveg_d13C) - 
      (tail(CO2conc_atmosphere,1) - CO2conc_atmosphere[1:par@yearsStage1])*par@CO2conc_fractionationFactor
    
    D13CSuessDifference_rhizodeposits[1:par@yearsStage1] <- 
      (tail(d13C_atmosphere,1) - vegData_stage1@d13C_rhizodeposits) - 
      (tail(CO2conc_atmosphere,1) - CO2conc_atmosphere[1:par@yearsStage1])*par@CO2conc_fractionationFactor
    
    # The d13C of vegetation in every year is calculated
    d13C_AGveg[1:par@yearsStage1] <- 
      d13C_atmosphere[1:par@yearsStage1] - D13CSuessDifference_AGveg[1:par@yearsStage1]
    d13C_BGveg[1:par@yearsStage1] <- 
      d13C_atmosphere[1:par@yearsStage1] - D13CSuessDifference_BGveg[1:par@yearsStage1]
    
    d13C_rhizodeposits[1:par@yearsStage1] <- 
      d13C_atmosphere[1:par@yearsStage1] - D13CSuessDifference_rhizodeposits[1:par@yearsStage1]
      
    if(par@nVeg == 2){
      
      # The difference in d13C between plant and d13CO2 is corrected for the [CO2],
      # based on Keeling et al. (2017)
      D13CSuessDifference_AGveg[(par@yearsStage1+1):par@numberOfSimulationYears] <- 
        (tail(d13C_atmosphere,1) - vegData_stage2@AGveg_d13C) - 
        (tail(CO2conc_atmosphere,1) - CO2conc_atmosphere[(par@yearsStage1+1):par@numberOfSimulationYears])*par@CO2conc_fractionationFactor
      D13CSuessDifference_BGveg[(par@yearsStage1+1):par@numberOfSimulationYears] <- 
        (tail(d13C_atmosphere,1) - vegData_stage2@BGveg_d13C) - 
        (tail(CO2conc_atmosphere,1) - CO2conc_atmosphere[(par@yearsStage1+1):par@numberOfSimulationYears])*par@CO2conc_fractionationFactor
      
      # The d13C of vegetation in every year is calculated
      d13C_AGveg[(par@yearsStage1+1):par@numberOfSimulationYears] <- 
        d13C_atmosphere[(par@yearsStage1+1):par@numberOfSimulationYears] - 
        D13CSuessDifference_AGveg[(par@yearsStage1+1):par@numberOfSimulationYears]
      d13C_BGveg[(par@yearsStage1+1):par@numberOfSimulationYears] <- 
        d13C_atmosphere[(par@yearsStage1+1):par@numberOfSimulationYears] - 
        D13CSuessDifference_BGveg[(par@yearsStage1+1):par@numberOfSimulationYears]
      
    } # End if(par@nVeg == 2)

  } else{ # If the Suess effect is not simulated
    
    # Write these codes...
    
  } # End if(par@Include_Suess_effect == 1)
  
  # The ratio of 13C/12C of atmospheric CO2
  atm_ratio_13C_12C <- ((d13C_atmosphere / 1000) + 1) * 0.0112372
  
  par@d13C_AGveg <- d13C_AGveg
  par@d13C_BGveg <- d13C_BGveg
  par@d13C_rhizodeposits <- d13C_rhizodeposits
  par@d13C_atmosphere <- d13C_atmosphere
  par@atm_ratio_13C_12C <- atm_ratio_13C_12C
  par@CO2conc_atmosphere <- CO2conc_atmosphere
    
  return(par)
    
}