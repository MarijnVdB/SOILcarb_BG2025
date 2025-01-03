# The annual inputs of 12C, 13C and 14C are constructed

constructIsotopeInputs <- function(par, vegData_stage1, vegData_stage2){
  
  # ----------------------
  # Aboveground vegetation
  # ----------------------
  
  # The 13C/12C ratio of vegetation over the simulated period
  ratio_13C_12C_atm <- ((par@d13C_atmosphere/1000)+1)*0.0112372
  AGveg_ratio_13C_12C <- ratio_13C_12C_atm*(1/((((par@d13C_atmosphere/1000)+1)/((par@d13C_AGveg/1000)+1))))
  
  # The factor that represents the fractionation against 13C by photosynthesis
  Fact <- ((par@d13C_AGveg/1000)+1)/((par@d13C_atmosphere/1000)+1)
  
  # The 14C/12C of leaves is calculated. The factor (Fact) is assumed
  # to be twice as large for 14C compared to 13C
  par@AGveg_ratio_14C_12C = par@atm_ratio_14C_12C*(1-(1-Fact)*2)
  
  # Annual inputs of 12C, 13C and 14C are calculated
  AGveg_12C_inputs <- rep(NA,par@numberOfSimulationYears)
  AGveg_13C_inputs <- rep(NA,par@numberOfSimulationYears)
  AGveg_14C_inputs <- rep(NA,par@numberOfSimulationYears)
  
  if(par@nVeg == 1){
    
    AGveg_12C_inputs <- vegData_stage1@i_agv / (1 + AGveg_ratio_13C_12C + par@AGveg_ratio_14C_12C)
    AGveg_13C_inputs <- (vegData_stage1@i_agv / AGveg_12C_inputs - par@AGveg_ratio_14C_12C-1) * AGveg_12C_inputs
    AGveg_14C_inputs <- vegData_stage1@i_agv - AGveg_12C_inputs - AGveg_13C_inputs
    
  } else if (par@nVeg == 2){
    
    AGveg_12C_inputs[1:par@yearsStage1] <- 
      vegData_stage1@i_agv / (1 + AGveg_ratio_13C_12C[1:yearsStage1] + par@AGveg_ratio_14C_12C[1:yearsStage1])
    AGveg_13C_inputs[1:par@yearsStage1] <- 
      (vegData_stage1@i_agv / AGveg_12C_inputs[1:yearsStage1] - par@AGveg_ratio_14C_12C[1:yearsStage1] - 1) * AGveg_12C_inputs[1:yearsStage1]
    AGveg_14C_inputs[1:par@yearsStage1] <- 
      vegData_stage1@i_agv - AGveg_12C_inputs[1:yearsStage1] - AGveg_13C_inputs[1:yearsStage1]
    
    AGveg_12C_inputs[(par@yearsStage1+1):par@numberOfSimulationYears] <- 
      vegData_stage2@i_agv / (1 + AGveg_ratio_13C_12C[(par@yearsStage1+1):par@numberOfSimulationYears] + 
      par@AGveg_ratio_14C_12C[(par@yearsStage1 + 1):par@numberOfSimulationYears])
    AGveg_13C_inputs[(par@yearsStage1+1):par@numberOfSimulationYears] <- 
      (vegData_stage2@i_agv / AGveg_12C_inputs[(par@yearsStage1 + 1):par@numberOfSimulationYears] - 
         par@AGveg_ratio_14C_12C[(par@yearsStage1 + 1):par@numberOfSimulationYears] - 1) * 
      AGveg_12C_inputs[(par@yearsStage1+1):par@numberOfSimulationYears]
    AGveg_14C_inputs[(par@yearsStage1+1):par@numberOfSimulationYears] <- 
      vegData_stage2@i_agv - AGveg_12C_inputs[(par@yearsStage1+1):par@numberOfSimulationYears] - 
      AGveg_13C_inputs[(par@yearsStage1+1):par@numberOfSimulationYears]
    
  } # End if(par@nVeg == 1)
  
  # ---------------------
  # Belowground vegetation
  # ---------------------
  
  # Inputs from dead roots
  
  BGveg_ratio_13C_12C_deadRoots <- ratio_13C_12C_atm*(1/((((par@d13C_atmosphere/1000)+1)/((par@d13C_BGveg/1000)+1))))
  
  # The factor that represents the fractionation against 13C by photosynthesis
  Fact <- (1/((((par@d13C_atmosphere/1000)+1)/((par@d13C_BGveg/1000)+1))))
  
  # The 14C/12C of roots is calculated, the factor (F) is assumed
  # to be twice as large for 14C compared to 13C
  BGveg_ratio_14C_12C_deadRoots <- par@atm_ratio_14C_12C*(1-(1-Fact)*2);
  
  # Inputs from rhizodeposition
  
  # The maximum d13C of AGveg
  max_d13C_AGveg <- max(par@d13C_AGveg)
  
  # An array to store the d13C of rhizodeposits in is created
  rhizodep_d13C_init <- rep(NA,par@numberOfSimulationYears)
  
  # If there are years with C4 vegetation, the d13C of rhizodeposits has to
  # be calculated differently for these years (no difference between d13C of 
  # aboveground vegetation and root exudates)
  
  # The years with C3 and C4 vegetation are determined
  r_C3 <- which(par@d13C_AGveg < -17)
  r_C4 <- which(par@d13C_AGveg >= -17)
  
  # The array with d13C of rhizodeposits is filled
  rhizodep_d13C_init[r_C3] <- par@d13C_rhizodeposits
  rhizodep_d13C_init[r_C4] <- par@d13C_AGveg[r_C4]       # Rhizodeposits have the same d13C as leaves
  
  
  BGveg_ratio_13C_12C_rhizodeposition <- ratio_13C_12C_atm*(1/((((par@d13C_atmosphere/1000)+1)/(((rhizodep_d13C_init)/1000)+1))))
  
  # The factor that represents the fractionation against 13C by photosynthesis
  Fact = (1/((((par@d13C_atmosphere/1000)+1)/(((rhizodep_d13C_init)/1000)+1))));
  
  # The 14C/12C of roots is calculated, the factor (F) is assumed
  # to be twice as large for 14C compared to 13C
  BGveg_ratio_14C_12C_rhizodeposition <- par@atm_ratio_14C_12C*(1-(1-Fact)*2);
  
  # Annual inputs of 12C, 13C and 14C are calculated
  BGveg_12C_inputs_deadRoots <- matrix(NA, nrow = par@layerNumber, ncol = par@numberOfSimulationYears)
  BGveg_13C_inputs_deadRoots <- matrix(NA, nrow = par@layerNumber, ncol = par@numberOfSimulationYears)
  BGveg_14C_inputs_deadRoots <- matrix(NA, nrow = par@layerNumber, ncol = par@numberOfSimulationYears)
  
  BGveg_12C_inputs_rhizodeposition <- matrix(NA, nrow = par@layerNumber, ncol = par@numberOfSimulationYears)
  BGveg_13C_inputs_rhizodeposition <- matrix(NA, nrow = par@layerNumber, ncol = par@numberOfSimulationYears)
  BGveg_14C_inputs_rhizodeposition <- matrix(NA, nrow = par@layerNumber, ncol = par@numberOfSimulationYears)
  
  # The dimensions of the relevant vectors are adjusted
  rootInputProfile <- vegData_stage1@i_bgveg_stage1
  dim(rootInputProfile) <- c(length(rootInputProfile),1)
  
  ratio_13C_12C_deadRoots <- BGveg_ratio_13C_12C_deadRoots
  dim(ratio_13C_12C_deadRoots) <- c(1,par@numberOfSimulationYears)
  
  ratio_14C_12C_deadRoots <- BGveg_ratio_14C_12C_deadRoots
  dim(ratio_14C_12C_deadRoots) <- c(1,par@numberOfSimulationYears)
  
  ratio_13C_12C_rhizodeposition <- BGveg_ratio_13C_12C_rhizodeposition
  dim(ratio_13C_12C_rhizodeposition) <- c(1,par@numberOfSimulationYears)
  
  ratio_14C_12C_rhizodeposition <- BGveg_ratio_14C_12C_rhizodeposition
  dim(ratio_14C_12C_rhizodeposition) <- c(1,par@numberOfSimulationYears)
  
  # Matrices of the correct dimensions for multiplication are constructed
  mat_rootInputProfile <- matrix(rep(rootInputProfile, each = par@yearsStage1), ncol = par@yearsStage1, byrow = TRUE)
  mat_ratio_13C_12C_deadRoots <- matrix(rep(ratio_13C_12C_deadRoots, each = par@layerNumber), nrow = par@layerNumber, byrow = FALSE)
  mat_ratio_14C_12C_deadRoots <- matrix(rep(ratio_14C_12C_deadRoots, each = par@layerNumber), nrow = par@layerNumber, byrow = FALSE)
  mat_ratio_13C_12C_rhizodeposition <- matrix(rep(ratio_13C_12C_rhizodeposition, each = par@layerNumber), nrow = par@layerNumber, byrow = FALSE)
  mat_ratio_14C_12C_rhizodeposition <- matrix(rep(ratio_14C_12C_rhizodeposition, each = par@layerNumber), nrow = par@layerNumber, byrow = FALSE)
  
  # Vegetation stage 1: inputs from dead roots
  BGveg_13C_inputs_deadRoots[,1:par@yearsStage1] <- (mat_rootInputProfile * (1 - par@Fdoc_rhizo)) / (1 + (1 / mat_ratio_13C_12C_deadRoots[,1:par@yearsStage1]))
  BGveg_12C_inputs_deadRoots[,1:par@yearsStage1] <- BGveg_13C_inputs_deadRoots[,1:par@yearsStage1] / mat_ratio_13C_12C_deadRoots[,1:par@yearsStage1]
  BGveg_14C_inputs_deadRoots[,1:par@yearsStage1] <- (mat_rootInputProfile * (1 - par@Fdoc_rhizo)) / (1 + 1 / (mat_ratio_14C_12C_deadRoots[,1:par@yearsStage1]))
  
  # Vegetation stage 1: inputs from rhizodeposits
  BGveg_13C_inputs_rhizodeposition[,1:par@yearsStage1] <- (mat_rootInputProfile * par@Fdoc_rhizo) / (1 + (1 / mat_ratio_13C_12C_rhizodeposition[,1:par@yearsStage1]))
  BGveg_12C_inputs_rhizodeposition[,1:par@yearsStage1] <- BGveg_13C_inputs_rhizodeposition[,1:par@yearsStage1] / mat_ratio_13C_12C_rhizodeposition[,1:par@yearsStage1]
  BGveg_14C_inputs_rhizodeposition[,1:par@yearsStage1] <- (mat_rootInputProfile * par@Fdoc_rhizo) / (1 + 1 / (mat_ratio_14C_12C_rhizodeposition[,1:par@yearsStage1]))
  
  # If there are 2 vegetation stages
  if(par@nVeg == 2){
    
    # The dimensions of the relevant vectors are adjusted
    rootInputProfile <- vegData_stage2@i_bgveg_stage2
    dim(rootInputProfile) <- c(length(rootInputProfile),1)
    
    # Matrices of the correct dimensions for multiplication are constructed
    mat_rootInputProfile <- matrix(rep(rootInputProfile, each = par@yearsStage2), ncol = par@yearsStage2, byrow = TRUE)
    
    # Vegetation stage 2: inputs from dead roots
    BGveg_13C_inputs_deadRoots[,(par@yearsStage1+1):par@numberOfSimulationYears] <- 
      (mat_rootInputProfile * (1 - par@Fdoc_rhizo)) / (1 + (1 / mat_ratio_13C_12C_deadRoots[,(par@yearsStage1+1):par@numberOfSimulationYears]))
    BGveg_12C_inputs_deadRoots[,(par@yearsStage1+1):par@numberOfSimulationYears] <- 
      BGveg_13C_inputs_deadRoots[,(par@yearsStage1+1):par@numberOfSimulationYears] /
      mat_ratio_13C_12C_deadRoots[,(par@yearsStage1+1):par@numberOfSimulationYears]
    BGveg_14C_inputs_deadRoots[,(par@yearsStage1+1):par@numberOfSimulationYears] <- 
      (mat_rootInputProfile * (1 - par@Fdoc_rhizo)) / (1 + (1 / mat_ratio_14C_12C_deadRoots[,(par@yearsStage1+1):par@numberOfSimulationYears]))
      
    # Vegetation stage 2: inputs from rhizodeposits
    BGveg_13C_inputs_rhizodeposition[,(par@yearsStage1+1):par@numberOfSimulationYears] <- 
      (mat_rootInputProfile * (1 - par@Fdoc_rhizo)) / (1 + (1 / mat_ratio_13C_12C_rhizodeposition[,(par@yearsStage1+1):par@numberOfSimulationYears]))
    BGveg_12C_inputs_rhizodeposition[,(par@yearsStage1+1):par@numberOfSimulationYears] <- 
      BGveg_13C_inputs_rhizodeposition[,(par@yearsStage1+1):par@numberOfSimulationYears] /
      mat_ratio_13C_12C_rhizodeposition[,(par@yearsStage1+1):par@numberOfSimulationYears]
    BGveg_14C_inputs_rhizodeposition[,(par@yearsStage1+1):par@numberOfSimulationYears] <- 
      (mat_rootInputProfile * (1 - par@Fdoc_rhizo)) / (1 + (1 / mat_ratio_14C_12C_rhizodeposition[,(par@yearsStage1+1):par@numberOfSimulationYears]))
  
  } # End if(par@nVeg == 2)
  
  # The ratio of 13C/12C is calculated for atmospheric CO2: this is necessary in the calculations
  ratioAtmosphic13CIsotopes <- ((par@d13C_atmosphere/1000)+1)*0.0112372;
  
  # Plotting
  # df <- data.frame("x" = 1:par@numberOfSimulationYears, "y" = BGveg_14C_inputs_rhizodeposition[1,])
  # ggplot(df, aes(x=x, y=y)) + 
  #   geom_line()
  
  # ------------------------------------------------------------------------------------------
  # The atmospheric CO2 concentration for 12CO2, 13CO2 and 14CO2 for every year is constructed
  # ------------------------------------------------------------------------------------------
  
  # The atmospheric CO2 concentration is converted to mol per cubic meter
  atmCO2Conc_MolPerCubicMeter<- (par@CO2conc_atmosphere * 42.2963) / 1e6
  
  atmCO2_12C_MolPerCubicMeter <- atmCO2Conc_MolPerCubicMeter / (1 + par@atm_ratio_13C_12C +par@atm_ratio_14C_12C)
  atmCO2_13C_MolPerCubicMeter <- (atmCO2Conc_MolPerCubicMeter / atmCO2_12C_MolPerCubicMeter - par@atm_ratio_14C_12C - 1) * atmCO2_12C_MolPerCubicMeter
  atmCO2_14C_MolPerCubicMeter <- atmCO2Conc_MolPerCubicMeter - atmCO2_12C_MolPerCubicMeter - atmCO2_13C_MolPerCubicMeter;
  
  # The isotopic inputs are combined in a list
  isotopicInputs <- list("AGveg_12C_inputs" = AGveg_12C_inputs,
                         "AGveg_13C_inputs" = AGveg_13C_inputs,
                         "AGveg_14C_inputs" = AGveg_14C_inputs,
                         "BGveg_12C_inputs_deadRoots" = BGveg_12C_inputs_deadRoots,
                         "BGveg_13C_inputs_deadRoots" = BGveg_13C_inputs_deadRoots,
                         "BGveg_14C_inputs_deadRoots" = BGveg_14C_inputs_deadRoots,
                         "BGveg_12C_inputs_rhizodeposition" = BGveg_12C_inputs_rhizodeposition,
                         "BGveg_13C_inputs_rhizodeposition" = BGveg_13C_inputs_rhizodeposition,
                         "BGveg_14C_inputs_rhizodeposition" = BGveg_14C_inputs_rhizodeposition)
  
  par@ratioAtmosphic13CIsotopes <- ratioAtmosphic13CIsotopes
  par@init_rhizo_13C_12C <- BGveg_ratio_13C_12C_rhizodeposition[1]
  par@init_rhizo_14C_12C <- BGveg_ratio_14C_12C_rhizodeposition[1]
  par@atmCO2_12C_MolPerCubicMeter <- atmCO2_12C_MolPerCubicMeter
  
  isotopeRatios <- list("AGveg_ratio_13C_12C" = AGveg_ratio_13C_12C,
                        "AGveg_ratio_14C_12C" = par@AGveg_ratio_14C_12C,
                        "BGveg_ratio_13C_12C_deadRoots" = BGveg_ratio_13C_12C_deadRoots,
                        "BGveg_ratio_14C_12C_deadRoots" = BGveg_ratio_14C_12C_deadRoots,
                        "BGveg_ratio_13C_12C_rhizodeposition" = BGveg_ratio_13C_12C_rhizodeposition,
                        "BGveg_ratio_14C_12C_rhizodeposition" = BGveg_ratio_14C_12C_rhizodeposition,
                        "BGveg_ratio_13C_12C_deadRoots" = BGveg_ratio_13C_12C_deadRoots,
                        "BGveg_ratio_14C_12C_deadRoots" = BGveg_ratio_14C_12C_deadRoots,
                        "BGveg_ratio_13C_12C_rhizodeposition" = BGveg_ratio_13C_12C_rhizodeposition,
                        "BGveg_ratio_14C_12C_rhizodeposition" = BGveg_ratio_14C_12C_rhizodeposition,
                        "atmCO2_12C_MolPerCubicMeter" = atmCO2_12C_MolPerCubicMeter,
                        "atmCO2_13C_MolPerCubicMeter" = atmCO2_13C_MolPerCubicMeter,
                        "atmCO2_14C_MolPerCubicMeter" = atmCO2_14C_MolPerCubicMeter)
  
  return(list("isotopicInputs" = isotopicInputs,
              "par" = par,
              "isotopeRatios" = isotopeRatios))
  
}