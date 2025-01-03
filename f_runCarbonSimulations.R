# This function initiates all variables before the SOILcarb model is run,
# runs the model and formats the necessary outputs

runCarbonSimulation <- function(calibParam, par1 = NA, vegData_stage1 = NA, vegData_stage2 = NA,
                                isotopicInputs = NA, isotopeRatios = NA, x_init = NA,
                                measuredData_litter = NA, measuredData_soil = NA){

  if (!is.na(calibParam[1])){ # If calibParam does not contain NAs, some data is grabbed from the global environment

    par1 = get("par", globalenv())
    vegData_stage1 = get("vegData_stage1", globalenv())
    vegData_stage2 = get("vegData_stage2", globalenv())
    isotopicInputs = get("isotopicInputs", globalenv())
    isotopeRatios = get("isotopeRatios", globalenv())
    measuredData_litter = get("measuredData_litter", globalenv())
    measuredData_soil = get("measuredData_soil", globalenv())
    x_init = get("x_init", globalenv())

  }

  # Some packages are loaded again (necessary in case of parallel calibration)
  if(par@calibMode == 1 | par@SAFE_sensitivity_isotopes == 1 | par@SAFE_sensitivity_allParameters == 1){
    library("deSolve")
    library("psych")

    source("f_various_functions.R")
    source("f_load_measured_data.R")
    source("f_inputParameters.R")
    source("f_load_d13C_data.R")
    source("f_load_d14C_data.R")
    source("f_construct_isotope_inputs.R")
    source("f_initialCO2DepthProfile.R")
    source("f_parameter_initialisation.R")
    source("f_runCarbonSimulations.R")
    source("f_SOILcarb.R")
    source("f_createClasses.R")
  }

  createInputClass() # Create the function to create the input object
  createVegClass()   # Create the function to create the vegetation object

  par <- par1

  initialStateVar <- x_init
 
  # ----------------------------------------------------------------------------
  # When calibrating, the calibration values are assigned to the correct 
  # parameters
  # ----------------------------------------------------------------------------
  
  if(par@calibMode == 1 | par@useCalibrationResults == 1 | par@runSeries == 1){

    if(par@site == "Hainich"){

      if(par@litterCalibration == 1){

        par@KmB_POC_L_stage1 <- calibParam[1]
        par@KmB_DOC_L_stage1 <- calibParam[2]
        par@Vmax_L <- calibParam[3]

      } else if (par@soilCalibration == 1) {
        
        par@VmaxD_root_R_stage1 <- calibParam[1]
        par@VmaxU_BioAv_stage1 <- calibParam[2]
        
        par@Db0_stage1 <- calibParam[3]
        
        vegData_stage1@betaRoots <- calibParam[4]
        vegData_stage1@fineRootBeta <- calibParam[4]
        
        par@Km_ads <- calibParam[5]
        par@Km_depol_M_stage1 <- calibParam[6]
        
        par@kDes_init <- calibParam[7]
        
        par@VmaxD_M_stage1 <- calibParam[8]
        par@Vmax_ads <- calibParam[9]

      } # End if(par@litterCalibration == 1)
    } # End if(par@site == "Hainich")

    # The parameters are initialized again
    out <- parameterInitialisation(par, vegData_stage1, vegData_stage2)
    par = out$par
    vegData_stage1 = out$vegData_stage1
    vegData_stage2 = out$vegData_stage2
    isotopicInputs = out$isotopicInputs
    isotopeRatios = out$isotopeRatios
    rm(out)

  } # End if(par@calibMode == 1) 
  else if (par@SAFE_sensitivity_isotopes == 1){ # In case the parameters for a SAFE sensitivity analysis have to be run
    
    vegData_stage1@AGveg_d13C <- as.numeric(calibParam[1])
    vegData_stage1@BGveg_d13C <- as.numeric(calibParam[2])
    vegData_stage1@d13C_rhizodeposits <- as.numeric(calibParam[3])
    par@fractIncorpCO2_soil <- as.numeric(calibParam[4])
    par@CO2conc_fractionationFactor <- as.numeric(calibParam[5])
    
    # The parameters are initialized again
    out <- parameterInitialisation(par, vegData_stage1, vegData_stage2)
    
    par = out$par
    vegData_stage1 = out$vegData_stage1
    vegData_stage2 = out$vegData_stage2
    isotopicInputs = out$isotopicInputs
    isotopeRatios = out$isotopeRatios
    
    rm(out)
    
  } else if (par@SAFE_sensitivity_allParameters == 1){  # In case the parameters for a SAFE sensitivity analysis have to be run
    
    par@Db0_stage1 <- as.numeric(calibParam[1])
    par@Db_eFold_depth_stage1 <- as.numeric(calibParam[2])
    vegData_stage1@betaRoots <- as.numeric(calibParam[3])
    vegData_stage1@fineRootBeta <- as.numeric(calibParam[3])
    
    par@VmaxD_root_R_stage1 <- as.numeric(calibParam[4])
    par@VmaxU_BioAv_stage1 <- as.numeric(calibParam[5])
    
    par@Vmax_ads <- as.numeric(calibParam[6])
    par@Km_ads <- as.numeric(calibParam[7])
    
    par@VmaxD_M_stage1 <- as.numeric(calibParam[8])
    par@Km_depol_M_stage1 <- as.numeric(calibParam[9])
    
    par@kDes_init <- as.numeric(calibParam[10])
    
    par@advectionRate_polyC_stage1 <- as.numeric(calibParam[11])
    
    # The parameters are initialized again
    out <- parameterInitialisation(par, vegData_stage1, vegData_stage2)
    
    par = out$par
    vegData_stage1 = out$vegData_stage1
    vegData_stage2 = out$vegData_stage2
    isotopicInputs = out$isotopicInputs
    isotopeRatios = out$isotopeRatios
    
    rm(out)
    
  }

  # ----------------------------------------------------------------------------
  # SOILcarb is run
  # ----------------------------------------------------------------------------
  
  currentStage <- 0  # A variable to keep track of the current vegetation stage

  for(ii in 1:par@nVeg){ # A loop over all vegetation stages
    
    currentStage <- currentStage + 1 # The currentStage variable is updated
    
    # Parameters values that depend on the vegetation stage are defined
    paramValues_currentStage <- stageParameters(par, vegData_stage1, 
                                                vegData_stage2, currentStage)
      
    # The initial state variables are formatted
    if(currentStage == 1){
      y_init <- initialStateVar
    } else if(currentStage == 2){
      # To do...
    }
    
    # The time sequence is defined
    if(currentStage == 1){
      times <- seq(1,par@yearsStage1)
    } else if(currentStage == 2){
      times <- seq(1,par@yearsStage2)
    }
    
    # The carbon inputs for the current vegetation stage are isolated
    if(currentStage == 1){
      c1 <- 1
      c2 <- par@yearsStage1
    } else if(currentStage == 2){
      c1 <- par@yearsStage1 + 1
      c2 <- par@yearsStage1 + par@yearsStage2
    }
    
    carbonInputs_currentStage <- list()   # An empty list
    
    # A loop to copy the relevant inputs for the years of this vegetation stage
    for(kk in 1:length(isotopicInputs)){
      if(dim(t(isotopicInputs[[kk]]))[1] == 1){   # If the matrix has 1 row
        carbonInputs_currentStage[[kk]] <- isotopicInputs[[kk]][c1:c2]
      } else if (nrow(isotopicInputs[[kk]]) > 1){ # If the matrix has >1 row
        carbonInputs_currentStage[[kk]] <- isotopicInputs[[kk]][,c1:c2]
      }
    }
    # The names are copied
    names(carbonInputs_currentStage) <- names(isotopicInputs)
    
    # The odes are solved
    out <- lsodes(func = SOILcarb, 
               y = y_init, 
               parms = paramValues_currentStage, 
               times = times, 
               carbonInputs = carbonInputs_currentStage,
               currentStage = currentStage,
               isotopeRatios = isotopeRatios,
               paramValues = par,
               isotopicInputs = isotopicInputs,
               hmax = 0)
    
    # The results are stored in variables
    if(currentStage == par@nVeg){
      outData <- formatOutput(out, 1950, vegData_stage1)
    }

    # --------------------------------------------------------------------------
    # If a calibration is run, the error is calculated
    # --------------------------------------------------------------------------
    
    # If a calibration is run, the errors are calculated. Also in case a series is run
    
    if(par@calibMode == 1 | par@runSeries == 1){
      
      res <- outData$res
      d13C <- outData$d13C
      d14C <- outData$d14C
      
      # Check if every year is simulated
      ndim <- length(res$litter_12C)

      # Check if no NA's are present in the output
      stop <-  0
      check <- F
      
      for(ii in 1:length(res)){
        check <- anyNA(res[ii], recursive = TRUE)
        if(check == T){
          stop <- 1
        }
      }
      
      for(ii in 1:length(d13C)){
        check <- anyNA(d13C[ii], recursive = TRUE)
        if(check == T){
          stop <- 1
        }
      }
      
      for(ii in 1:length(d14C)){
        check <- anyNA(d14C[ii], recursive = TRUE)
        if(check == T){
          stop <- 1
        }
      }
      
      if(ndim < par@numberOfSimulationYears | stop == 1){ # If not all necessary years have been run
        
        relError_tot <- 1e3
        
      } else {
        
      # ------------------------------
      # The errors for the soil layers
      # ------------------------------
      
      # The row numbers for which modelled values have to be averaged are retrieved
      lowerRows <- (measuredData_soil$lowerBound * 100) + 1
      upperRows <- measuredData_soil$upperBound * 100
      
      # The modelled C depth profiles are converted to %OC
      POC_rhizo_tot <- (res$POC_R[dim(res$POC_R)[1],] / par@soilMassArray) * 100
      
      # All pools in the bulk soil are assumed to be mineral-associated C
      bulSoilCarbon <- res$MIN_B[dim(res$MIN_B)[1],] + res$POLY_B[dim(res$POLY_B)[1],] + res$MIC_B[dim(res$MIC_B)[1],]
      bulkSoilCarbon_tot  <-  (bulSoilCarbon / par@soilMassArray) * 100
      Cmic_rhizo_tot <- (res$MIC_R[dim(res$MIC_R)[1],] / par@soilMassArray) * 100
      Cmic_bulk_tot <- (res$MIC_B[dim(res$MIC_B)[1],] / par@soilMassArray) * 100
      
      # The d13C and d14C of the bulk soil C are calculated
      endOf <- dim(res$MIN_B)[1]
      
      bulkSoilCarbon_d13C <- ((d13C$MIN_B[endOf,] *res$MIN_B[endOf,]) + (d13C$POLY_B[endOf,] * res$POLY_B[endOf,]) + (d13C$MIC_B[endOf,] *res$MIC_B[endOf,])) / bulSoilCarbon
      bulkSoilCarbon_d14C <- ((d14C$MIN_B[endOf,] *res$MIN_B[endOf,]) + (d14C$POLY_B[endOf,] * res$POLY_B[endOf,]) + (d14C$MIC_B[endOf,] *res$MIC_B[endOf,])) / bulSoilCarbon
      
      # The model results are interpolated to match the measurement depths
      
      # The boundary depths are changed to multiples of 0.01
      boundaryDepth_tmp <- floor(par@boundaryDepth * 100) / 100
      maxInterpDepth <- floor(tail(par@midDepth_layer,1) * 100) / 100
      
      # Model results are interpolated to 0.01 m depth intervals
      interpDepths <- seq(0.01, maxInterpDepth, 0.01)
      endOf <- dim(res$MIN_B)[1]
      POC_rhizo_tot_interp <- approx(par@midDepth_layer,POC_rhizo_tot,interpDepths)$y
      bulkSoilCarbon_tot_interp <- approx(par@midDepth_layer,bulkSoilCarbon_tot,interpDepths)$y
      Cmic_rhizo_tot_interp <- approx(par@midDepth_layer,Cmic_rhizo_tot,interpDepths)$y
      Cmic_bulk_tot_interp <- approx(par@midDepth_layer,Cmic_bulk_tot,interpDepths)$y
      POC_d13C_interp <- approx(par@midDepth_layer,d13C$POC_R[endOf,],interpDepths)$y
      POC_d14C_interp <- approx(par@midDepth_layer,d14C$POC_R[endOf,],interpDepths)$y
      bulkSoilCarbon_d13C_interp <- approx(par@midDepth_layer,bulkSoilCarbon_d13C,interpDepths)$y
      bulkSoilCarbon_d14C_interp <- approx(par@midDepth_layer,bulkSoilCarbon_d14C,interpDepths)$y
      d14C_soil_interp <- approx(par@midDepth_layer,d14C$soilC[endOf,],interpDepths)$y
      
      # These are averaged for the measurement rows so they can be compared
      # to the measurements, together with the other modeled values
      POC_OC_averaged <- rep(NA,length(lowerRows))
      bulkSoilCarbon_OC_averaged <- rep(NA,length(lowerRows))
      Cmic_rhizo_averaged <- rep(NA,length(lowerRows))
      Cmic_bulk_averaged <- rep(NA,length(lowerRows))
      POC_d13C_averaged <- rep(NA,length(lowerRows))
      POC_d14C_averaged <- rep(NA,length(lowerRows))
      bulkSoilCarbon_d13C_averaged <- rep(NA,length(lowerRows))
      bulkSoilCarbon_d14C_averaged <- rep(NA,length(lowerRows))
      totalSOC_14C_averaged <- rep(NA,length(lowerRows))
      
      # The modeled C at the measurement depth are calculated
      for(i in 1:length(lowerRows)){
        
        POC_OC_averaged[i] <- mean(POC_rhizo_tot_interp[lowerRows[i]:upperRows[i]])
        bulkSoilCarbon_OC_averaged[i] <- mean(bulkSoilCarbon_tot_interp[lowerRows[i]:upperRows[i]])
        Cmic_rhizo_averaged[i] <- mean(Cmic_rhizo_tot_interp[lowerRows[i]:upperRows[i]])
        Cmic_bulk_averaged[i] <- mean(Cmic_bulk_tot_interp[lowerRows[i]:upperRows[i]])
        POC_d13C_averaged[i] <- mean(POC_d13C_interp[lowerRows[i]:upperRows[i]])
        POC_d14C_averaged[i] <- mean(POC_d14C_interp[lowerRows[i]:upperRows[i]])
        bulkSoilCarbon_d13C_averaged[i] <- mean(bulkSoilCarbon_d13C_interp[lowerRows[i]:upperRows[i]])
        bulkSoilCarbon_d14C_averaged[i] <- mean(bulkSoilCarbon_d14C_interp[lowerRows[i]:upperRows[i]])
        totalSOC_14C_averaged[i] <- mean(d14C_soil_interp[lowerRows[i]:upperRows[i]])
        
      }
      
      # ---------------------------------------------
      # Normalized sum of squares for the soil layers
      # ---------------------------------------------
      
      # Soil
      # POC
      n <- length(POC_OC_averaged)
      range_C = max(measuredData_soil$fPOC_conc[1:n]) - min(measuredData_soil$fPOC_conc[1:n])
      relError_OC_POC <- sum(((measuredData_soil$fPOC_conc[1:n] - POC_OC_averaged) / range_C)^2)

      # MAOC
      n <- length(bulkSoilCarbon_OC_averaged)
      range_C = max(measuredData_soil$MAOC_conc[1:n]) - min(measuredData_soil$MAOC_conc[1:n])
      relError_OC_MAOC <- sum(((measuredData_soil$MAOC_conc[1:n] - bulkSoilCarbon_OC_averaged) / range_C)^2)
      
      # d13C POC
      n <- length(POC_d13C_averaged)
      range_C = max(measuredData_soil$fPOC_d13C[1:n]) - min(measuredData_soil$fPOC_d13C[1:n])
      relError_d13C_POC <- sum(((measuredData_soil$fPOC_d13C[1:n] - POC_d13C_averaged) / range_C)^2)
      
      # d13C MAOC
      n <- length(bulkSoilCarbon_d13C_averaged)
      range_C = max(measuredData_soil$MAOC_d13C[1:n]) - min(measuredData_soil$MAOC_d13C[1:n])
      relError_d13C_MAOC <- sum(((measuredData_soil$MAOC_d13C[1:n] - bulkSoilCarbon_d13C_averaged) / range_C)^2)
      
      # d14C POC
      NA_rows <- !is.na(measuredData_soil$fPOC_d14C) # TRUE for the rows containing data, FALSE for rows containing NA
      range_C = max(measuredData_soil$fPOC_d14C[NA_rows]) - min(measuredData_soil$fPOC_d14C[NA_rows])
      relError_d14C_POC <- sum(((measuredData_soil$fPOC_d14C[NA_rows] - POC_d14C_averaged[NA_rows]) / range_C)^2)
      
      # d14C MAOC
      NA_rows <- !is.na(measuredData_soil$MAOC_d14C) # TRUE for the rows containing data, FALSE for rows containing NA
      range_C = max(measuredData_soil$MAOC_d14C[NA_rows]) - min(measuredData_soil$MAOC_d14C[NA_rows])
      relError_d14C_MAOC <- sum(((measuredData_soil$MAOC_d14C[NA_rows] - bulkSoilCarbon_d14C_averaged[NA_rows]) / range_C)^2)
      
      # ---------------------------
      # Errors for the litter layer
      # ---------------------------
      
      endOf <- length(res$C_litter)
      relError_litter_OC <- ((measuredData_litter$Litter_Cstock - res$C_litter[endOf]) / measuredData_litter$Litter_Cstock)^2
      relError_litter_DOC <- ((measuredData_litter$Litter_Cstock*(0.33) - res$DOC_L[endOf]) / measuredData_litter$Litter_Cstock*(0.33))^2
      relError_litter_POC <- ((measuredData_litter$Litter_Cstock*(0.66) - res$POC_L[endOf]) / measuredData_litter$Litter_Cstock*(0.66))^2
      relError_litter_MIC <- ((measuredData_litter$Litter_Cstock*(0.01) - res$MIC_L[endOf]) / measuredData_litter$Litter_Cstock*(0.01))^2
      relError_litter_d13C <- ((measuredData_litter$Litter_d13C - d13C$C_L[endOf]) / measuredData_litter$Litter_d13C)^2
      relError_litter_d14C <- ((measuredData_litter$Litter_d14C - d14C$C_L[endOf]) / measuredData_litter$Litter_d14C)^2
      
      } # End if(ndim < par@numberOfSimulationYears)
 
    } # End calculating errors
    
    # The final errors are calculated, based on the calibration scenario
    if(par@calibMode == 1 | par@runSeries == 1){
      
      if(par@site == "Hainich"){

        # Some parameter combinations cause the solver to stop prematurely. In this case, the error is given a large value
        # The dimensions of the results
        ndim <- length(res$litter_12C)
        if(ndim < par@numberOfSimulationYears | stop == 1){ # If not all necessary years have been run
          relError_tot <- 1e3
        } else { # If all necessary years have been run
        
        # ---------------------------
        # Errors for the litter layer
        # ---------------------------
        
        if (par@calibCheck$C_d13C_d14C == 1){
        
          # --- C, d13C, d14C ---
          relError_litter <- sum(c(relError_litter_OC, relError_litter_DOC, relError_litter_POC, relError_litter_d13C, relError_litter_d14C))
        
        } else if (par@calibCheck$C_d13C == 1){
        
          # --- C, d13C ---
          relError_litter <- sum(c(relError_litter_OC, relError_litter_DOC, relError_litter_POC, relError_litter_d13C))
          
        } else if (par@calibCheck$C == 1) {
        
          # --- C ---
          relError_litter <- sum(c(relError_litter_OC, relError_litter_DOC, relError_litter_POC))
        
        } else if (par@calibCheck$C_d14C == 1){
          
          # --- C, d14C ---
          relError_litter <- sum(c(relError_litter_OC, relError_litter_DOC, relError_litter_POC, relError_litter_d14C))
        
        }
          
        # -------------------
        # Errors for the soil
        # -------------------
        
        if (par@calibCheck$C_d13C_d14C == 1){
        
          # --- C, d13C, d14C ---
          relError_soil <- sum(c(relError_OC_POC, relError_OC_MAOC, relError_d13C_POC,
                                 relError_d13C_MAOC, relError_d14C_POC, relError_d14C_MAOC))
        
        } else if (par@calibCheck$C_d13C == 1){
        
          # --- C, d13C ---
          relError_soil <- sum(c(relError_OC_POC, relError_OC_MAOC, relError_d13C_POC, relError_d13C_MAOC))
        
        } else if (par@calibCheck$C == 1){
        
          # --- C ---
          relError_soil <- sum(c(relError_OC_POC, relError_OC_MAOC))
        
        } else if (par@calibCheck$C_d14C == 1){
        
          # --- C, d14C ---
          relError_soil <- sum(c(relError_OC_POC, relError_OC_MAOC, relError_d14C_MAOC, relError_d14C_POC))
        
        }
        
        # -----------------------------
        # The final errors are assigned
        # -----------------------------
        
        if(par@litterCalibration == 1){
          relError_tot <- relError_litter
        } else if (par@soilCalibration == 1){
          relError_tot <- relError_soil
        }
        
        # ------------------------
        # Constraints for the soil
        # ------------------------
        
        if(par@soilCalibration == 1){
          
          # The mass of BIOAV should be smaller than POC
          bioAv_mass <- sum(res$BIOAV_R[endOf,])
          POC_mass <- sum(res$POC_R[endOf,])
          if(bioAv_mass > POC_mass){
            relError_tot <- 1e3
          }
          
          # The mass of soil microbes cannot be >5% of total OC
          Ctot_soil_mass <- sum(res$C_soil[endOf,])
          Cmic_soil_mass <- sum(res$MIC_R[endOf,] + res$MIC_B[endOf,])
          if(Cmic_soil_mass > Ctot_soil_mass*0.05){
            relError_tot <- 1e3
          }
          
          # The turnover time of rhizosphere and bulk soil microbes
          # should be within certain ranges
          
          rowNum <- dim(res[["MIC_R"]])[1]
          Ctot_R <- res[["MIC_R"]][rowNum,] + res[["BIOAV_R"]][rowNum,] + res[["POC_R"]][rowNum,]
          K_R <- par@K_mic_R * Ctot_R
          MicToBioAv_R <- res[["MIC_R"]][rowNum,] / res[["BIOAV_R"]][rowNum,]
          uptake_R <- par@VmaxU_BioAv_stage1 * res[["BIOAV_R"]][rowNum,] * (MicToBioAv_R / (par@KmU_Bioav_profile_stage1 + MicToBioAv_R))
          growthRate_R <- (uptake_R * par@CUE_R) / res[["MIC_R"]][rowNum,]
          death_R <- (growthRate_R * res[["MIC_R"]][rowNum,]^2) / K_R # Abs. amount of death microbes per year
          TT_rhizo_yr <- res[["MIC_R"]][rowNum,] / death_R # Tunrover time = stock / flux [year]
          TT_rhizo_day <- TT_rhizo_yr * 365 # Turnover time converted to days
          
        } # End if(par@soilCalibration == 1)
        
      } # End if(par@site == "Hainich")
      } # End if(ndim < par@numberOfSimulationYears)
    }  # End if(par@calibMode == 1)
    
    # --------------------------------------------------------------
    # In case the parameters for a SAFE sensitivity analysis are run
    # --------------------------------------------------------------
    
    if(par@SAFE_sensitivity_isotopes == 1){
      
      d13C <- outData$d13C
      d14C <- outData$d14C
      
      output <- rep(NA,6)
      
      # 1. d13C at the soil surface
      output[1] <- d13C$soilC[nrow(d13C$soilC),1]
      # 2. d13C in the 8th layer (midDepth = 0.41 cm)
      output[2] <- d13C$soilC[nrow(d13C$soilC),8]
      # 3. d14C at the soil surface
      output[3] <- d14C$soilC[nrow(d14C$soilC),1]
      # 4. d14C in the 8th layer (midDepth = 0.41 cm)
      output[4] <- d14C$soilC[nrow(d14C$soilC),8]
      # 5. The difference in d13C between the surface and ca. 40 cm depth
      output[5] <- output[2] - output[1]
      # 6. The difference in d14C between the surface and ca. 40 cm depth
      output[6] <- output[4] - output[3]

    }
    
    if(par@SAFE_sensitivity_allParameters == 1){
      
      res <- outData$res
      d13C <- outData$d13C
      d14C <- outData$d14C
      
      output <- rep(NA,7)
      
      # 1. d13C at the soil surface
      output[1] <- d13C$soilC[nrow(d13C$soilC),1]
      # 2. d13C in the 8th layer (midDepth = 0.41 cm)
      output[2] <- d13C$soilC[nrow(d13C$soilC),8]
      # 3. d14C at the soil surface
      output[3] <- d14C$soilC[nrow(d14C$soilC),1]
      # 4. d14C in the 8th layer (midDepth = 0.41 cm)
      output[4] <- d14C$soilC[nrow(d14C$soilC),8]
      # 5. The difference in d13C between the surface and ca. 40 cm depth
      output[5] <- output[2] - output[1]
      # 6. The difference in d14C between the surface and ca. 40 cm depth
      output[6] <- output[4] - output[3]
      
      # 7. Carbon stock
      Cstock <- sum(res$C_soil[par@numberOfSimulationYears,])
      output[7] <- Cstock
    }
    
    # ----------------------------------------
    # The output of the function is determined
    # ----------------------------------------

    if(par@calibMode == 1){
      return(relError_tot)
    } else if (par@runSeries == 1){
      return(c(outData, "error" = relError_tot,
               "relError_OC_POC" = relError_OC_POC,
               "relError_OC_MAOC" = relError_OC_MAOC,
               "relError_d13C_POC" = relError_d13C_POC,
               "relError_d13C_MAOC" = relError_d13C_MAOC,
               "relError_d14C_POC" = relError_d14C_POC,
               "relError_d14C_MAOC" = relError_d14C_MAOC))
    } else if(par@SAFE_sensitivity_isotopes == 1){
      return(output)
    } else if (par@SAFE_sensitivity_allParameters == 1){
      return(output)
    }else if (par@normalRun == 1){
      return(outData)
    }
    
  }
  
}