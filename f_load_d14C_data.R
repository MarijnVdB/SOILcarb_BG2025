# This script is used to load the 14C input data

loadD14CData <- function(par){
  
  # The data is loaded
  allData_d14C = read.csv('d14CData.csv')
  
  # Only data for the simulated period is needed
  startyear_14C <- par@lastRealSimulationYear - par@numberOfSimulationYears + 1;
  endyear_14C <- par@lastRealSimulationYear;
  
  # Only data for the simulated period is needed
  startyear_14C <- par@lastRealSimulationYear - par@numberOfSimulationYears + 1
  endyear_14C <- par@lastRealSimulationYear
  
  # Determine the start- and endrow of the d14CO2 data that's needed
  if(par@numberOfSimulationYears > nrow(allData_d14C)){
    # If more years are simulated than d14CO2 years are available
    startRow <- 1
    endRow <- which(allData_d14C$Year == endyear_14C)
  } else{ # If less years are simulated than d14CO2 years are available
    startRow <- which(allData_d14C$Year == startyear_14C)
    endRow <- which(allData_d14C$Year == endyear_14C)
  }
  
  # The time series are adjusted accordingly
  Atmospheric_d14C_data <- rep(NA,par@numberOfSimulationYears)
  Atmospheric_d14C_years <- rep(NA,par@numberOfSimulationYears)
  
  if(par@numberOfSimulationYears > nrow(allData_d14C)){
    # If the number of simulated years is larger than the number of years
    # for which atmospheric 14C data is available
    firstProvidedYear_atmD14C <- nrow(allData_d14C)
    
    endOf <- length(Atmospheric_d14C_data)
    Atmospheric_d14C_data[(endOf - (firstProvidedYear_atmD14C) - 1):endOf] <- allData_d14C$d14CO2
    Atmospheric_d14C_data[1:(endOf - firstProvidedYear_atmD14C)] <- allData_d14C$d14CO2[1]
    
    Atmospheric_d14C_years <- seq(par@lastRealSimulationYear - par@numberOfSimulationYears + 1,par@lastRealSimulationYear)
    
  } else { # If the number of simulated years is less than the number of years with d14C data for atmospheric CO2
    Atmospheric_d14C_data <- allData_d14C$d14CO2[startRow:endRow]
    Atmospheric_d14C_years <- allData_d14C$Year[startRow:endRow]
  }
  
  # The 14C/12C ratio of atmospheric CO2 is calculated, for a d13C of -25 per mil
  atm_ratio_14C_12C_corr = ((Atmospheric_d14C_data/1000)+1)*(0.95*(1.18*10^(-12))*exp((1950-1950)/8267))
  # And corrected for the d13C value of atmospheric CO2
  atm_ratio_14C_12C = atm_ratio_14C_12C_corr/(0.9750/(1+((par@d13C_atmosphere)/1000)))^2
  
  par@atm_ratio_14C_12C <- atm_ratio_14C_12C
  par@d14C_atmosphere <- Atmospheric_d14C_data
  par@d14C_atmosphere_years <- Atmospheric_d14C_years
  
  return(par)
  
} # End function