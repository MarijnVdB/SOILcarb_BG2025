# The core SOILcarb model

SOILcarb <- function(Time, State, Pars, carbonInputs, currentStage, isotopeRatios, paramValues, isotopicInputs){
  # Time is the time
  # State contains the state variables
  # Pars contains the parameter values for the current stage
  # carbonInputs contains the annual carbon inputs
  # p contains the general model parameters
  
  if(Time > paramValues@numberOfSimulationYears){
    Time  <-  paramValues@numberOfSimulationYears
  }
  
  yr <- floor(Time)   # The current year
  
  # The number of simulation years left is displayed
  if((paramValues@numberOfSimulationYears - yr) %% 100 == 0 & paramValues@calibMode == 0 & paramValues@runSeries == 0 & paramValues@SAFE_sensitivity_isotopes == 0 & paramValues@SAFE_sensitivity_allParameters == 0){
    print(paste("Still", as.character(paramValues@numberOfSimulationYears - yr), "years to go"))
  }
  
  # The d13C and d14C value of atmospheric CO2 is retrieved
  atm_d13C <- paramValues@d13C_atmosphere[yr]
  atm_d14C <- paramValues@d14C_atmosphere[yr]
  
  # The inputs of 12C, 13C and 14C for the current year are retrieved
  i_12L <- carbonInputs$AGveg_12C_inputs[yr]
  i_13L <- carbonInputs$AGveg_13C_inputs[yr]
  i_14L <- carbonInputs$AGveg_14C_inputs[yr]
  
  i_12R_deathRoots <- carbonInputs$BGveg_12C_inputs_deadRoots[,yr]
  i_13R_deathRoots <- carbonInputs$BGveg_13C_inputs_deadRoots[,yr]
  i_14R_deathRoots <- carbonInputs$BGveg_14C_inputs_deadRoots[,yr]

  i_12R_rhizoDepo <-carbonInputs$BGveg_12C_inputs_rhizodeposition[,yr]
  i_13R_rhizoDepo <- carbonInputs$BGveg_13C_inputs_rhizodeposition[,yr]
  i_14R_rhizoDepo <- carbonInputs$BGveg_14C_inputs_rhizodeposition[,yr]
  
  # The concentrations of 12CO2, 13CO2 and 14CO2 of atmospheric CO2 are retrieved
  allYears <- seq(1,paramValues@numberOfSimulationYears)
  if(currentStage == 1){
    k <- yr
  } else if (currentStage == 2){
    k <- paramValues@yearsStage1 + yr
  }
  
  i_atm12CO2 <- isotopeRatios$atmCO2_12C_MolPerCubicMeter[k]
  i_atm13CO2 <- isotopeRatios$atmCO2_13C_MolPerCubicMeter[k]
  i_atm14CO2 <- isotopeRatios$atmCO2_14C_MolPerCubicMeter[k]
  rm(k, allYears)
  
  # --------------------------------------------------------
  # The state variables are assigned to increase readability
  # --------------------------------------------------------
  
  v <- list(
  
  # === 12C ===
  "POC_12L" = State["POC_12L"],        # Litter POC12
  "DOC_12L" = State["DOC_12L"],       # Litter DOC12
  "MIC_12L" = State["MIC_12L"],       # Litter MIC12

  # === 13C ===
  "POC_13L" = State["POC_13L"],        # Litter POC13
  "DOC_13L" = State["DOC_13L"],      # Litter DOC13
  "MIC_13L" = State["MIC_13L"],       # Litter MIC13

  # === 14C ===
  "POC_14L" = State["POC_14L"],       # Litter POC13
  "DOC_14L" = State["DOC_14L"],       # Litter DOC13
  "MIC_14L" = State["MIC_14L"],       # Litter MIC13

  # Soil - Rhizosphere
  
  # === 12C ===
  "POC_12R"  = State[startsWith(names(State), "POC_12R")],    # Rhizosphere POC12
  "BIOAV_12R"  = State[startsWith(names(State), "BIOAV_12R")],  # Rhizosphere BIOAV12
  "MIC_12R"  = State[startsWith(names(State), "MIC_12R")],     # Rhizosphere MIC12
  "CO2_12R"  = State[startsWith(names(State), "CO2_12R")],     # Rhizosphere 12CO2
  
  # === 13C ===
  "POC_13R"  = State[startsWith(names(State), "POC_13R")],     # Rhizosphere POC13
  "BIOAV_13R"  = State[startsWith(names(State), "BIOAV_13R")],  # Rhizosphere BIOAV13
  "MIC_13R"  = State[startsWith(names(State), "MIC_13R")],     # Rhizosphere MIC13
  "CO2_13R"  = State[startsWith(names(State), "CO2_13R")],     # Rhizosphere 13CO2
  
  # === 14C ===
  "POC_14R"  = State[startsWith(names(State), "POC_14R")],     # Rhizosphere POC14
  "BIOAV_14R"  = State[startsWith(names(State), "BIOAV_14R")],  # Rhizosphere BIOAV14
  "MIC_14R"  = State[startsWith(names(State), "MIC_14R")],    # Rhizosphere MIC14
  "CO2_14R"  = State[startsWith(names(State), "CO2_14R")],     # Rhizosphere 14CO2
  
  # Soil - Bulk soil
  
  # === 12C ===
  "POLY_12B"  = State[startsWith(names(State), "POLY_12B")],    # Bulk soil POLY_12
  "MIN_12B"  =  State[startsWith(names(State), "MIN_12B")],    # Bulk soil MIN_12
  "MIC_12B"  = State[startsWith(names(State), "MIC_12B")],     # Bulk soil MIC12
  "CO2_12B"  = State[startsWith(names(State), "CO2_12B")],     # Bulk soil 12CO2
  "SURF_12B"  = State[startsWith(names(State), "SURF_12B")],     # Bulk soil SURF_12
  
  # === 13C ===
  "POLY_13B"  = State[startsWith(names(State), "POLY_13B")],    # Bulk soil POLY_13
  "MIN_13B"  = State[startsWith(names(State), "MIN_13B")],     # Bulk soil MIN_13
  "MIC_13B"  = State[startsWith(names(State), "MIC_13B")],    # Bulk soil MIC13
  "CO2_13B"  = State[startsWith(names(State), "CO2_13B")],    # Bulk soil 13CO2
  "SURF_13B"  = State[startsWith(names(State), "SURF_13B")],     # Bulk soil SURF_12
  
  # === 14C ===
  "POLY_14B"  = State[startsWith(names(State), "POLY_14B")],    # Bulk soil POLY_14
  "MIN_14B"  = State[startsWith(names(State), "MIN_14B")],    # Bulk soil MIN_14
  "MIC_14B"  = State[startsWith(names(State), "MIC_14B")],     # Bulk soil MIC14
  "CO2_14B"  = State[startsWith(names(State), "CO2_14B")],     # Bulk soil 14CO2
  "SURF_14B"  = State[startsWith(names(State), "SURF_14B")]     # Bulk soil SURF_12
  
  ) # Close list
  
  # --------------------------------------------
  # The d13C and d14C of soil CO2 are calculated
  # --------------------------------------------
  
  CO2_12 = v$CO2_12R + v$CO2_12B
  CO2_13 = v$CO2_13R + v$CO2_13B
  CO2_14 = v$CO2_14R + v$CO2_14B
  
  d13C_CO2 <- calculateD13C(CO2_12,CO2_13)
  d14C_CO2 <- calculateD14C(CO2_12,CO2_14,d13C_CO2,1950)
  
  # -----------------------
  # The litter layer is run
  # -----------------------
  
  # === 12C ===
  
  # Inputs to POC
  F1_12 <- i_12L * (1 - Pars$L_fractionLeachable)
  
  # Inputs to DOC
  F2_12 <- i_12L * Pars$L_fractionLeachable
  
  # Depolimerisation of POC
  F3_12 <- (Pars$Vmax_L * v$POC_12L * v$MIC_12L) / (Pars$KmB_POC_L * (1 + (v$POC_12L/Pars$KmB_POC_L) + (v$DOC_12L/Pars$KmB_DOC_L)) + v$MIC_12L)
  
  # Depolimerisation of DOC
  F4_12 <- (Pars$Vmax_L * v$DOC_12L * v$MIC_12L) / (Pars$KmB_DOC_L * (1 + (v$POC_12L/Pars$KmB_POC_L) + (v$DOC_12L/Pars$KmB_DOC_L)) + v$MIC_12L)
  
  # Death of microbes
  Ctot_L <- v$MIC_12L + v$DOC_12L + v$POC_12L
  K_L <- Pars$K_mic_L * Ctot_L
  growthRate_L <- (F3_12 * Pars$CUE_L + F4_12 * Pars$CUE_L) / v$MIC_12L
  death_L <- (growthRate_L * v$MIC_12L^2) / K_L
  
  # Death microbial biomass to POC
  F5_12 <- death_L * (1-Pars$solubleMicFraction)
  
  # Death microbial biomass to DOC
  F6_12 <- death_L * Pars$solubleMicFraction
  
  # Flux of C from POC and DOC to microbes
  F7_12 <- F3_12 + F4_12
  
  # Microbial uptake of 12C, 13C and 14C from atmospheric CO2
  atmUptake_12C_L <- F7_12 * Pars$CUE_L * Pars$fractIncorpCO2_litter
  atmUptake_13C_L <-  ((atm_d13C / 1000) + 1) * 0.0112372 * atmUptake_12C_L
  atmCO2_ratio_14C_12C_25 <- ((atm_d14C / 1000) + 1) * (0.95 * (1.18 * 10^(-12)) * exp((1950-1950)/8267))
  atmCO2_ratio_14C_12C <- atmCO2_ratio_14C_12C_25 / ((0.9750 / (1 + (atm_d13C/1000)))^2)
  atmUptake_14C_L <- atmCO2_ratio_14C_12C * atmUptake_12C_L
  
  # The fluxes from C to microbes are reduced according to the amount of CO2-C that is taken up
  F3_12 <- F3_12 - (F3_12 * Pars$CUE_L * Pars$fractIncorpCO2_litter)
  F4_12 <- F4_12 - (F4_12 * Pars$CUE_L * Pars$fractIncorpCO2_litter)
  
  # C losses from the litter layer through bioturbation and leaching
  F8_12 <- (v$POC_12L - F3_12) * Pars$bioturbFract # C lost through bioturbation
  CO2_loss <- F7_12 * (1 - Pars$CUE_L) #Total C lost as CO2
  F9_12 <- (v$DOC_12L - F4_12) * Pars$leachFract # C lost through leaching
  
  # === 13C ===
  
  # Inputs to POC
  F1_13 <- i_13L * (1 - Pars$L_fractionLeachable)
  
  # Inputs to DOC
  F2_13 <- i_13L * Pars$L_fractionLeachable
  
  # Depolimerisation of POC
  F3_13 <- v$POC_13L * (F3_12 / v$POC_12L)
  
  # Depolimerisation of DOC
  F4_13 <- v$DOC_13L * (F4_12 / v$DOC_12L)
  
  # Death microbial biomass to POC
  F5_13 <- v$MIC_13L * (F5_12 / v$MIC_12L)
  
  # Death microbial biomass to DOC
  F6_13 <- v$MIC_13L * (F6_12 / v$MIC_12L)
  
  # Flux of C from POC and DOC to microbes, including heterotrophic CO2 assimilation
  F7_13 <- F3_13 + F4_13
  
  # C losses from the litter layer through bioturbation and leaching
  F8_13 <- (v$POC_13L - F3_13) * Pars$bioturbFract
  F9_13 <- v$DOC_13L * (F9_12 / v$DOC_12L)
  
  # === 14C ===
  
  # Inputs to POC
  F1_14 <- i_14L * (1 - Pars$L_fractionLeachable)
  
  # Inputs to DOC
  F2_14 <- i_14L * Pars$L_fractionLeachable
  
  # Depolimerisation of POC
  F3_14 <- v$POC_14L * (F3_12 / v$POC_12L)
  
  # Depolimerisation of DOC
  F4_14 <- v$DOC_14L * (F4_12 / v$DOC_12L)
  
  # Death microbial biomass to POC
  F5_14 <- v$MIC_14L * (F5_12 / v$MIC_12L)
  
  # Death microbial biomass to DOC
  F6_14 <- v$MIC_14L * (F6_12 / v$MIC_12L)
  
  # Flux of C from POC and DOC to microbes, including heterotrophic CO2 assimilation
  F7_14 <- F3_14 + F4_14
  
  # C losses from the litter layer through bioturbation and leaching
  F8_14 <- (v$POC_14L - F3_14) * Pars$bioturbFract
  F9_14 <- v$DOC_14L * (F9_12 / v$DOC_12L)
  
  # ----------------------------
  # The soil is run: rhizosphere
  # ----------------------------

  # === 12C ===
  
  # Inputs from death roots to the rhizosphere POC pool
  F10_12 <- i_12R_deathRoots
  
  # Inputs of litter POC to the rhizosphere POC pool
  Fpoc_12 <- rep(0, Pars$layerNumber)
  Fpoc_12[1] <- F8_12
  
  # Inputs from rhizodeposition in the rhizosphere bio-available C pool
  F11_12 <- i_12R_rhizoDepo
  
  # Depolimerisation of rhizosphere POC to bio-available C
  MicToPoc_R <- v$MIC_12R /  v$POC_12R
  F12_12 <- Pars$VmaxD_root_R * v$POC_12R * (MicToPoc_R / (Pars$KmB_root_R_profile + MicToPoc_R))
  
  # Uptake of bio-available C by rhizosphere microbes (not accounting for CUE)
  MicToBioAv_R <- v$MIC_12R / v$BIOAV_12R
  F13_12 <- Pars$VmaxU_BioAv * v$BIOAV_12R * (MicToBioAv_R / (Pars$KmU_Bioav_profile + MicToBioAv_R))
  
  # An error is generated when the amount of bio-availabe C taken up by
  # microbes is larger then the amount available
  if(any(F13_12 > v$BIOAV_12R)){
    print("Warning: the amount of BIOAV taken up is larger than the amount available")
    # F13_12 <-  v$BIOAV_12R * 0.8
  }
  
  # Transfer of bio-available C to polymeric C
  FbioToPoly_12 <- Pars$bioToPoly * (v$BIOAV_12R - F13_12)
  
  # Microbial death in the rhizosphere
  Ctot_R <- v$MIC_12R + v$BIOAV_12R + v$POC_12R
  K_R <- Pars$K_mic_R * Ctot_R
  # K_R <- Pars$K_mic_R * v$BIOAV_12R
  growthRate_R <- (F13_12 * Pars$CUE_R) / v$MIC_12R
  F14_12 <- (growthRate_R * v$MIC_12R^2) / K_R
  
  # Microbial uptake of 12C, 13C and 14C from soil CO2
  atmUptake_12C_R <- F13_12 * Pars$CUE_R * Pars$fractIncorpCO2_soil
  atmUptake_13C_R <- (((d13C_CO2 / 1000) + 1) * 0.0112372) * atmUptake_12C_R
  atmCO2_ratio_14C_12C_25 <- (((d14C_CO2 / 1000) + 1) * (0.95 * (1.18 * 10^(-12)) * exp((1950 - 1950) / 8267)))
  atmCO2_ratio_14C_12C <- atmCO2_ratio_14C_12C_25 / ((0.9750 / (1 + (d13C_CO2 / 1000)))^2)
  atmUptake_14C_R <- atmCO2_ratio_14C_12C * atmUptake_12C_R
  
  # === 13C ===
  
  # Inputs from death roots to the rhizosphere POC pool
  F10_13 <- i_13R_deathRoots
  
  # Inputs of litter POC to the rhizosphere POC pool
  Fpoc_13 <- rep(0, Pars$layerNumber)
  Fpoc_13[1] <- F8_13
  
  # Inputs from rhizodeposition in the rhizosphere bio-available C pool
  F11_13 <- i_13R_rhizoDepo
  
  # Depolimerisation of rhizosphere POC to bio-available C
  F12_13 <- v$POC_13R * (F12_12 / v$POC_12R)
  
  # Uptake of bio-available C by rhizosphere microbes
  F13_13 <- v$BIOAV_13R * (F13_12 / v$BIOAV_12R)

  # Transfer of bio-available C to polymeric C
  FbioToPoly_13 <- v$BIOAV_13R * (FbioToPoly_12 / v$BIOAV_12R)
  
  # Microbial death in the rhizosphere
  F14_13 <- v$MIC_13R * (F14_12 / v$MIC_12R)
  
  # === 14C ===
  
  # Inputs from death roots to the rhizosphere POC pool
  F10_14 <- i_14R_deathRoots
  
  # Inputs of litter POC to the rhizosphere POC pool
  Fpoc_14 <- rep(0, Pars$layerNumber)
  Fpoc_14[1] <- F8_14
  
  # Inputs from and rhizodeposition in the rhizosphere bio-available C pool
  F11_14 <- i_14R_rhizoDepo
  
  # Depolimerisation of rhizosphere POC to bio-available C
  F12_14 <- v$POC_14R * (F12_12 / v$POC_12R)
  
  # Uptake of bio-availble C by rhizosphere microbes
  F13_14 <- v$BIOAV_14R * (F13_12 / v$BIOAV_12R)

  # Transfer of bio-available C to polymeric C
  FbioToPoly_14 <- v$BIOAV_14R * (FbioToPoly_12 / v$BIOAV_12R)
  
  # Microbial death in the rhizosphere
  F14_14 <- v$MIC_14R * (F14_12 / v$MIC_12R)
  
  # --------------------------
  # The soil is run: bulk soil
  # --------------------------
  
  # === 12C ===
  
  # Protection of polyC in the bulk soil
  surf <- v$SURF_12B * Pars$rhizosphereVolumeFraction
  F15_12 <- (Pars$Vmax_ads * v$POLY_12B * surf) / ((Pars$Km_ads * (1 + (surf/Pars$Km_ads) + (v$MIC_12B/Pars$Km_depol_M))) + v$POLY_12B)
  
  # De-protection of MAOC in the bulk soil
  F16_12 <- Pars$kDes * v$MIN_12B
  
  # Depolimerisation and uptake of polymeric C by microbes in the bulk soil
  F17_12 <- (Pars$VmaxD_M * v$POLY_12B * v$MIC_12B) / ((Pars$Km_depol_M * (1 + (surf/Pars$Km_ads) + (v$MIC_12B/Pars$Km_depol_M))) + v$POLY_12B)
  
  # Microbial death in the bulk soil
  Ctot_B <- v$MIC_12B + v$POLY_12B + v$MIN_12B
  K_B <- Pars$K_mic_B * Ctot_B
  growthRate_B <- (F17_12 * Pars$CUE_B) / v$MIC_12B
  F18_12 <- (growthRate_B * v$MIC_12B^2) / K_B
  
  # Inputs of litter DOC to the polyC pool
  F19_12 <- rep(0,Pars$layerNumber)
  F19_12[1] <- F9_12
  
  # Microbial uptake of 12C, 13C and 14C from soil CO2
  atmUptake_12C_B <- F17_12 * Pars$CUE_B * Pars$fractIncorpCO2_soil
  atmUptake_13C_B <- (((d13C_CO2 / 1000) + 1) * 0.0112372) * atmUptake_12C_B
  atmCO2_ratio_14C_12C_25 <- (((d14C_CO2 / 1000) + 1) * (0.95 * (1.18 * 10^(-12)) * exp((1950-1950)/8267)))
  atmCO2_ratio_14C_12C <- atmCO2_ratio_14C_12C_25 / ((0.9750 / (1 + (d13C_CO2/1000)))^2)
  atmUptake_14C_B <- atmCO2_ratio_14C_12C * atmUptake_12C_B
  
  # The fluxes from C to microbes is reduced according to the amount of CO2-C is taken up
  F17_12 <- F17_12 - (F17_12 * Pars$CUE_B * Pars$fractIncorpCO2_soil)
  
  # === 13C ===
  
  # Protection of polyC in the bulk soil
  F15_13 <- v$POLY_13B * (F15_12 / v$POLY_12B)
  
  # De-protection of MAOC in the bulk soil
  F16_13 <- v$MIN_13B * (F16_12 / v$MIN_12B)
  
  # Depolimerisation and uptake of polymeric C by microbes in the bulk soil
  F17_13 <- v$POLY_13B * (F17_12 / v$POLY_12B)
  
  # Microbial death in the bulk soil
  F18_13 <- v$MIC_13B * (F18_12 / v$MIC_12B)
  
  # Inputs of litter DOC to the polyC pool
  F19_13 <- rep(0,Pars$layerNumber)
  F19_13[1] <- F9_13
  
  # === 14C ===
  
  # Protection of polyC in the bulk soil
  F15_14 <- v$POLY_14B * (F15_12 / v$POLY_12B)
  
  # De-protection of MAOC in the bulk soil
  F16_14 <- v$MIN_14B * (F16_12 / v$MIN_12B)
  
  # Depolimerisation and uptake of polymeric C by microbes in the bulk soil
  F17_14 <- v$POLY_14B * (F17_12 / v$POLY_12B)
  
  # Microbial death in the bulk soil
  F18_14 <- v$MIC_14B * (F18_12 / v$MIC_12B)
  
  # Inputs of litter DOC to the polyC pool
  F19_14 <- rep(0,Pars$layerNumber)
  F19_14[1] <- F9_14
  
  # ------------------------------------
  # Bioturbation
  # ------------------------------------
  # https://mathformeremortals.wordpress.com/2013/01/12/a-numerical-second-derivative-from-three-points/
  
  depth_tmp <- c(0, Pars$midDepth_layer, tail(Pars$boundaryDepth,1))
  
  # If bioturbation is simulated
  if(Pars$includeBioturbation == 1){
    
    # ----------------------------
    # Methods Goffin et al. (2014)
    # ----------------------------
    
    Db <- c(Pars$DbArray[1], Pars$DbArray, tail(Pars$DbArray,1))
    endOf <- length(paramValues@midDepth_layer)
    depth_bioturb <- c(paramValues@midDepth_layer[1] - (paramValues@midDepth_layer[2] - paramValues@midDepth_layer[1]),
                       paramValues@midDepth_layer,
                       paramValues@midDepth_layer[endOf] + (paramValues@midDepth_layer[endOf] - paramValues@midDepth_layer[endOf-1]) 
                       )

    # The profiles containing carbon data
    allC <- matrix(c(v$POC_12R,
                     v$POC_13R,
                     v$POC_14R,
                     
                     v$POLY_12B,
                     v$POLY_13B,
                     v$POLY_14B,
                     
                     v$MIN_12B,
                     v$MIN_13B,
                     v$MIN_14B,
                     
                     v$MIC_12B,
                     v$MIC_13B,
                     v$MIC_14B,
                     
                     F19_12,
                     F19_13,
                     F19_14,
                     
                     Fpoc_12,
                     Fpoc_13,
                     Fpoc_14,
                      
                     v$SURF_12B,
                     v$SURF_13B,
                     v$SURF_14B),
                   nrow = 21, byrow = TRUE)
    
    # The deltaC matrix is calculated
    endOf <- dim(allC)[2]
    allC_withBounds <- matrix(c(allC[,1], allC, allC[,endOf]), nrow = 21, byrow = FALSE)
    
    # The thickness of each layer
    endOf <- length(Pars$boundaryDepth)
    dz_layer <- Pars$boundaryDepth[2:endOf] - Pars$boundaryDepth[1:(endOf-1)]
    dz_layer <- c(dz_layer[1], dz_layer, tail(dz_layer,1))
    
    # The C concentration of each layer (since it's a diffusive process)
    s <- dim(allC_withBounds)[1]
    allC_withBounds <- allC_withBounds / matrix(rep(dz_layer,s), nrow = s, byrow = TRUE)
    
    # Ftop
    endOf <- length(Db)
    D <- harmonic.mean(rbind(Db[1:(endOf-2)], Db[2:(endOf-1)]))
    endOf <- dim(allC_withBounds)[2]
    C1 <- allC_withBounds[,1:(endOf-2)]
    C2 <- allC_withBounds[,2:(endOf-1)]
    endOf <- length(depth_bioturb)
    dz <- depth_bioturb[2:(endOf-1)] - depth_bioturb[1:(endOf-2)]
    s <- dim(allC_withBounds)[1]
    Ftop <- -matrix(rep(D,s), nrow = s, byrow = TRUE) * (C1-C2) / matrix(rep(dz,s), nrow = s, byrow = TRUE)
    
    # Fbot
    endOf <- length(Db)
    D <- harmonic.mean(rbind(Db[2:(endOf-1)], Db[3:endOf]))
    endOf <- dim(allC_withBounds)[2]
    C1 <- allC_withBounds[,2:(endOf-1)]
    C2 <- allC_withBounds[,3:endOf]
    endOf <- length(depth_bioturb)
    dz <- depth_bioturb[3:endOf] - depth_bioturb[2:(endOf-1)]
    Fbot <- -matrix(rep(D,s), nrow = s, byrow = TRUE) * (C1-C2) / matrix(rep(dz,s), nrow = s, byrow = TRUE)
    
    # The change in C concentration in each layer
    endOf <- length(Pars$boundaryDepth)
    dz <- Pars$boundaryDepth[2:endOf] - Pars$boundaryDepth[1:(endOf-1)]
    dC <- -(Ftop - Fbot) / matrix(rep(dz,s), nrow = s, byrow = TRUE)
    
    # This is converted to the change in total C
    endOf <- length(dz_layer)
    dC <- dC * matrix(rep(dz_layer[2:(endOf-1)],s), nrow = s, byrow = TRUE)
    
    # The fluxes per C pool and isotope are isolated
    Fb1 <- dC[1,]        # POC_12R
    Fb2 <- dC[2,]        # POC_13R
    Fb3 <- dC[3,]        # POC_14R
    Fb4 <- dC[4,]        # POLY_12B
    Fb5 <- dC[5,]        # POLY_13B
    Fb6 <- dC[6,]        # POLY_14B
    Fb7 <- dC[7,]        # MIN_12B
    Fb8 <- dC[8,]        # MIN_13B
    Fb9 <- dC[9,]        # MIN_14B
    Fb10 <- dC[10,]      # MIC_12B
    Fb11 <- dC[11,]      # MIC_13B
    Fb12 <- dC[12,]      # MIC_14B
    Fb13 <- dC[13,]      # Leached DOC_12
    Fb14 <- dC[14,]      # Leached DOC_13
    Fb15 <- dC[15,]      # Leached DOC_14
    Fb16 <- dC[16,]      # Bioturbated POC_12
    Fb17 <- dC[17,]      # Bioturbated POC_13
    Fb18 <- dC[18,]      # Bioturbated POC_14
    Fb19 <- dC[19,]      # SURF_12B
    Fb20 <- dC[20,]      # SURF_13B
    Fb21 <- dC[21,]      # SURF_14B
    
  } else { # If bioturbation is not simulated
    
    Fb1 <- rep(0,Pars$layerNumber)        # POC_12R
    Fb2 <- rep(0,Pars$layerNumber)        # POC_13R
    Fb3 <- rep(0,Pars$layerNumber)        # POC_14R
    Fb4 <- rep(0,Pars$layerNumber)        # POLY_12B
    Fb5 <- rep(0,Pars$layerNumber)        # POLY_13B
    Fb6 <- rep(0,Pars$layerNumber)        # POLY_14B
    Fb7 <- rep(0,Pars$layerNumber)        # MIN_12B
    Fb8 <- rep(0,Pars$layerNumber)        # MIN_13B
    Fb9 <- rep(0,Pars$layerNumber)        # MIN_14B
    Fb10 <- rep(0,Pars$layerNumber)       # MIC_12B
    Fb11 <- rep(0,Pars$layerNumber)       # MIC_13B
    Fb12 <- rep(0,Pars$layerNumber)       # MIC_14B
    Fb13 <- rep(0,Pars$layerNumber)       # Leached DOC_12
    Fb14 <- rep(0,Pars$layerNumber)       # Leached DOC_13
    Fb15 <- rep(0,Pars$layerNumber)       # Leached DOC_14
    Fb16 <- rep(0,Pars$layerNumber)       # Bioturbated POC_12
    Fb17 <- rep(0,Pars$layerNumber)       # Bioturbated POC_12
    Fb18 <- rep(0,Pars$layerNumber)       # Bioturbated POC_12
    Fb19 <- rep(0,Pars$layerNumber)       # SURF_12B
    Fb20 <- rep(0,Pars$layerNumber)       # SURF_13B
    Fb21 <- rep(0,Pars$layerNumber)       # SURF_14B
    
  }

  # ------------------------------------
  # Advection
  # -----------------------------------

  if(Pars$includeAdvection){
    
    # The profiles containing carbon data
    allC <- matrix(c(v$BIOAV_12R,
                     v$BIOAV_13R,
                     v$BIOAV_14R,
                     
                     v$POLY_12B,
                     v$POLY_13B,
                     v$POLY_14B),
                   
                   nrow = 6, byrow = TRUE)
    
    s <- dim(allC)[1]
    allC_withBounds <- matrix(c(rep(0,s), allC), nrow = 6, byrow = FALSE)
    
    # The thickness of each layer
    endOf <- length(Pars$boundaryDepth)
    dz_layer <- Pars$boundaryDepth[2:endOf] - Pars$boundaryDepth[1:(endOf-1)]
    dz_layer <- c(dz_layer[1], dz_layer)
    
    # The C concentration of each layer
    allC_withBounds <- allC_withBounds / matrix(rep(dz_layer,s), nrow = s, byrow = TRUE)
    
    # The amount of C coming in from the top of the layer
    endOf <- dim(allC_withBounds)[2]
    Ftop <- Pars$advectionRate_polyC * allC_withBounds[,1:(endOf-1)]
    
    # The amount of C leaving the bottom of the layer
    Fbot <- Pars$advectionRate_polyC * allC_withBounds[,2:endOf]
    
    # The change in concentration
    Fall <- Ftop - Fbot
    
    Fa1 <- Fall[1,] # BIOAV_12R
    Fa2 <- Fall[2,] # BIOAV_13R
    Fa3 <- Fall[3,] # BIOAV_14R
    
    Fa4 <- Fall[4,] # POLY_12B
    Fa5 <- Fall[5,] # POLY_13B
    Fa6 <- Fall[6,] # POLY_14B
    
  } else { # If no advection is simulated
    
    Fa1 <- rep(0,Pars$layerNumber) # BIOAV_12R
    Fa2 <- rep(0,Pars$layerNumber) # BIOAV_13R
    Fa3 <- rep(0,Pars$layerNumber) # BIOAV_14R
    
    Fa4 <- rep(0,Pars$layerNumber) # POLY_12B
    Fa5 <- rep(0,Pars$layerNumber) # POLY_13B
    Fa6 <- rep(0,Pars$layerNumber) # POLY_14B
    
  }
  
  # ----------------
  # Diffusion of CO2
  # ----------------
  
  # The C losses from decomposition are converted from kg C to mole CO2 per cubic meter
  C_losses_kgCperyear_R <- unname(rbind(F13_12 * (1 - Pars$CUE_R), F13_13 * (1 - Pars$CUE_R), F13_14 * (1 - Pars$CUE_R)))
  C_losses_kgCperyear_B <- unname(rbind(F17_12 * (1 - Pars$CUE_B), F17_13 * (1 - Pars$CUE_B), F17_14 * (1 - Pars$CUE_B)))
  
  C_Losses_gramCO2PerDay_B <- ((C_losses_kgCperyear_R * 1e3) / 12) * 44   # gram CO2 per layer
  C_Losses_gramCO2PerDay_R <- ((C_losses_kgCperyear_B * 1e3) / 12) * 44   # gram CO2 per layer
  
  C_Losses_MolePerDay_B <- C_Losses_gramCO2PerDay_B / 44   # Mole CO2 per soil layer
  C_Losses_MolePerDay_R <- C_Losses_gramCO2PerDay_R / 44   # Mole CO2 per soil layer
  
  C_Losses_MolePerCubicMeterPerDay_R <- C_Losses_MolePerDay_R * (1 / matrix(rep(t(Pars$layerThickness),3), nrow = 3, byrow = TRUE)) * (1 / matrix(rep(t(Pars$airFilledPorosityDepthProfile),3), nrow = 3, byrow = TRUE))   # Mole per cubic meter
  C_Losses_MolePerCubicMeterPerDay_B <- C_Losses_MolePerDay_B *( 1 / matrix(rep(t(Pars$layerThickness),3), nrow = 3, byrow = TRUE)) * (1 / matrix(rep(t(Pars$airFilledPorosityDepthProfile),3), nrow = 3, byrow = TRUE))   # Mole per cubic meter
  
  # The annual CO2 production profile is converted from gram C to mole CO2 per cubic meter
  annualRootCO2Production_gramCO2 <- ((Pars$annualRootCO2 * 1e3) / 12) * 44   # gram CO2 per layer
  annualRootCO2Production_Mole <- annualRootCO2Production_gramCO2 / 44   # Mole CO2 per soil layer
  annualRootCO2Production_MolePerCubicMeter <- annualRootCO2Production_Mole * (1 / Pars$layerThickness) * (1 / Pars$airFilledPorosityDepthProfile)   # Mole per cubic meter
  
  # The 12C, 13C and 14C of root inputs for the evaluated year are calculated
  d13C_roots <- calculateD13C(isotopicInputs$BGveg_12C_inputs_deadRoots[1,yr], isotopicInputs$BGveg_13C_inputs_deadRoots[1,yr])  # The d13C of roots in the evaluated year
  ratio_14C_12C_roots <- isotopicInputs$BGveg_14C_inputs_deadRoots[1,yr] / isotopicInputs$BGveg_12C_inputs_deadRoots[1,yr]
  
  anualRootCO2Production_13C_12C <- ((d13C_roots / 1000) + 1) * 0.0112372
  annualRootCO2Production_14C_12C <- ratio_14C_12C_roots
  annualRootCO2Production_12C <- annualRootCO2Production_MolePerCubicMeter / (1 + anualRootCO2Production_13C_12C + annualRootCO2Production_14C_12C)
  annualRootCO2Production_13C <- (annualRootCO2Production_MolePerCubicMeter / annualRootCO2Production_12C - annualRootCO2Production_14C_12C - 1) * annualRootCO2Production_12C
  annualRootCO2Production_14C <- annualRootCO2Production_MolePerCubicMeter - annualRootCO2Production_12C - annualRootCO2Production_13C
  
  # Inputs of CO2, units are mole per cubic meter
  P_12C_R <- annualRootCO2Production_12C + C_Losses_MolePerCubicMeterPerDay_R[1,]
  P_13C_R <- annualRootCO2Production_13C + C_Losses_MolePerCubicMeterPerDay_R[2,]
  P_14C_R <- annualRootCO2Production_14C + C_Losses_MolePerCubicMeterPerDay_R[3,]
  
  P_12C_B <- C_Losses_MolePerCubicMeterPerDay_B[1,]
  P_13C_B <- C_Losses_MolePerCubicMeterPerDay_B[2,]
  P_14C_B <- C_Losses_MolePerCubicMeterPerDay_B[3,]
  
  conc_12R <- unname(c(i_atm12CO2, v$CO2_12R, tail(v$CO2_12R,1)))
  conc_13R <- unname(c(i_atm13CO2, v$CO2_13R, tail(v$CO2_13R,1)))
  conc_14R <- unname(c(i_atm14CO2, v$CO2_14R, tail(v$CO2_14R,1)))
  
  conc_12B <- unname(c(i_atm12CO2, v$CO2_12B, tail(v$CO2_12B,1)))
  conc_13B <- unname(c(i_atm13CO2, v$CO2_13B, tail(v$CO2_13B,1)))
  conc_14B <- unname(c(i_atm14CO2, v$CO2_14B, tail(v$CO2_14B,1)))
  
  CO2conc_R = matrix(c(conc_12R, conc_13R, conc_14R), nrow = 3, byrow = TRUE)
  CO2conc_B = matrix(c(conc_12B, conc_13B, conc_14B), nrow = 3, byrow = TRUE)
  
  endOf <- length(depth_tmp)
  x1 <- depth_tmp[1:(endOf-2)]
  x2 <- depth_tmp[2:(endOf-1)]
  x3 <- depth_tmp[3:endOf]
  rm(endOf)
  
  endOf <- dim(CO2conc_R)[2]
  y1R <- CO2conc_R[,1:(endOf-2)]
  y2R <- CO2conc_R[,2:(endOf-1)]
  y3R <- CO2conc_R[,3:endOf]
  rm(endOf)
  
  endOf <- dim(CO2conc_B)[2]
  y1B <- CO2conc_B[,1:(endOf-2)]
  y2B <- CO2conc_B[,2:(endOf-1)]
  y3B <- CO2conc_B[,3:endOf]
  rm(endOf)
  
  d2_12R <- ((2 * y1R[1,]) / ((x2 - x1) * (x3 - x1))) - ((2 * y2R[1,]) / ((x3 - x2) * (x2 - x1))) + ((2 * y3R[1,]) / ((x3 - x2) * (x3-x1)))
  d2_13R <- ((2 * y1R[2,]) / ((x2 - x1) * (x3 - x1))) - ((2 * y2R[2,]) / ((x3 - x2) * (x2 - x1))) + ((2 * y3R[2,]) / ((x3 - x2) * (x3-x1)))
  d2_14R <- ((2 * y1R[3,]) / ((x2 - x1) * (x3 - x1))) - ((2 * y2R[3,]) / ((x3 - x2) * (x2 - x1))) + ((2 * y3R[3,]) / ((x3 - x2) * (x3-x1)))
  
  d2_12B <- ((2 * y1B[1,]) / ((x2 - x1) * (x3 - x1))) - ((2 * y2B[1,]) / ((x3 - x2) * (x2 - x1))) + ((2 * y3B[1,]) / ((x3-x2) * (x3-x1)))
  d2_13B <- ((2 * y1B[2,]) / ((x2 - x1) * (x3 - x1))) - ((2 * y2B[2,]) / ((x3 - x2) * (x2 - x1))) + ((2 * y3B[2,]) / ((x3-x2) * (x3-x1)))
  d2_14B <- ((2 * y1B[3,]) / ((x2 - x1) * (x3 - x1))) - ((2 * y2B[3,]) / ((x3 - x2) * (x2 - x1))) + ((2 * y3B[3,]) / ((x3-x2) * (x3-x1)))
  
  Fc_12R <- paramValues@De_12C_year * d2_12R + P_12C_R
  Fc_13R <- paramValues@De_13C_year * d2_13R + P_13C_R
  Fc_14R <- paramValues@De_14C_year * d2_14R + P_14C_R
  
  Fc_12B <- paramValues@De_12C_year * d2_12B + P_12C_B
  Fc_13B <- paramValues@De_13C_year * d2_13B + P_13C_B
  Fc_14B <- paramValues@De_14C_year * d2_14B + P_14C_B
  
  # -------------------------------------
  # The dydt output vector is constructed
  # -------------------------------------
  
  # For radioactive decay of 14C
  dt <- 1     # Time step is 1 year 
  fractionLeft <- exp(-(1.21e-4)*dt)
  
  dydt <- list(c(
    
    # Litter 12C
    unname(F1_12 + F5_12 - F3_12 - F8_12), # POC_12L
    unname(F2_12 + F6_12 - F4_12 - F9_12), # DOC_12L
    unname(F3_12 * Pars$CUE_L + F4_12 * Pars$CUE_L + atmUptake_12C_L - F5_12 - F6_12), # MIC_12L
    
    # Litter 13C
    unname(F1_13 + F5_13 - F3_13 - F8_13), # POC_13L
    unname(F2_13 + F6_13 - F4_13 - F9_13), # DOC_13L
    unname(F3_13 * Pars$CUE_L + F4_13 * Pars$CUE_L + atmUptake_13C_L - F5_13 - F6_13), # MIC_13L
    
    # Litter 14C
    unname(F1_14 + F5_14 - F3_14 - F8_14 - v$POC_14L * (1 - fractionLeft)), # POC_14L
    unname(F2_14 + F6_14 - F4_14 - F9_14 - v$DOC_14L *(1 - fractionLeft)), # DOC_14L
    unname(F3_14 * Pars$CUE_L + F4_14 * Pars$CUE_L + atmUptake_14C_L - F5_14 - F6_14 - v$MIC_14L * (1 - fractionLeft)), # MIC_14L
    
    # Rhizosphere 12C
    unname(F10_12 - F12_12 + Fpoc_12 + Fb1 + Fb16), # POC_12R
    unname(F11_12 + F12_12 + F14_12 * Pars$solubleMicFraction - (F13_12 - atmUptake_12C_R) - FbioToPoly_12 + Fa1), # BIOAV_12R 
    unname((F13_12 - atmUptake_12C_R) * Pars$CUE_R + atmUptake_12C_R - F14_12), # MIC_12R
    unname(Fc_12R), # CO2_12R
    
    # Rhizosphere 13C
    unname(F10_13 - F12_13 + Fpoc_13 + Fb17 + Fb2), # POC_13R 
    unname(F11_13 + F12_13 + F14_13 * Pars$solubleMicFraction - (F13_13 - atmUptake_13C_R) - FbioToPoly_13 + Fa2), # BIOAV_13R 
    unname((F13_13 - atmUptake_13C_R) * Pars$CUE_R + atmUptake_13C_R - F14_13), # MIC_13R
    unname(Fc_13R), # CO2_13R
    
    # Rhizosphere 14C
    unname(F10_14 - F12_14 + Fpoc_14 - v$POC_14R * (1 - fractionLeft) + Fb3 + Fb18), # POC_14R 
    unname(F11_14 + F12_14 + F14_14 * Pars$solubleMicFraction - (F13_14 - atmUptake_14C_R) - FbioToPoly_14 - v$BIOAV_14R * (1 - fractionLeft) + Fa3), # BIOAV_14R 
    unname((F13_14 - atmUptake_14C_R) * Pars$CUE_R + atmUptake_14C_R - F14_14 - v$MIC_14R * (1 - fractionLeft)), # MIC_14R
    unname(Fc_14R), # CO2_14R
    
    # Bulk soil 12C
    unname(F14_12 * (1 - Pars$solubleMicFraction) + F16_12 +
             F18_12 + F19_12 - F15_12 - F17_12 + Fb4 + Fb13 + Fa4 + FbioToPoly_12), # POLY_12B 
    unname(F15_12 - F16_12 + Fb7), # MIN_12B 
    unname(F17_12 * Pars$CUE_B + atmUptake_12C_B - F18_12 + Fb10), # MIC_12B  
    unname(Fc_12B), # CO2_12B
    unname(F16_12 - F15_12 + Fb19), # SURF_12B
    
    # Bulk soil 13C
    unname(F14_13 * (1 - Pars$solubleMicFraction) + F16_13 +
             F18_13 + F19_13 - F15_13 - F17_13 + Fb5 + Fb14 + Fa5 + FbioToPoly_13), # POLY_13B 
    unname(F15_13 - F16_13 + Fb8), # MIN_13B 
    unname(F17_13 * Pars$CUE_B + atmUptake_13C_B - F18_13 + Fb11), # MIC_13B  
    unname(Fc_13B), # CO2_13B
    unname(F16_13 - F15_13 + Fb20), # SURF_12B
    
    # Bulk soil 14C
    unname(F14_14 * (1 - Pars$solubleMicFraction) + F16_14 +
             F18_14 + F19_14 - F15_14 - F17_14 - v$POLY_14B * (1 - fractionLeft) + Fb6 + Fb15  + Fa6 + FbioToPoly_14), # POLY_14B 
    unname(F15_14 - F16_14 - v$MIN_14B * (1 - fractionLeft) + Fb9), # MIN_14B 
    unname(F17_14 * Pars$CUE_B + atmUptake_14C_B - F18_14 -
             v$MIC_14B * (1 - fractionLeft) + Fb12), # MIC_14B 
    unname(Fc_14B),  # CO2_14B
    unname(F16_14 - F15_14 + Fb21) # SURF_12B
    
  )
  )
  
  return(dydt)
  
}