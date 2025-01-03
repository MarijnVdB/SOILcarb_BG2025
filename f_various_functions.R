# ------------------------------------------------------------------------------
# Calculate the thickness of the simulated layers
# ------------------------------------------------------------------------------

calculateLayers <-function(layer1_thickness, maxDepth, boundaryFactor){

  layerThickness  <-  layer1_thickness;
  boundaryDepth  <-  layerThickness;
  midDepth  <-  layerThickness/2;

  i = 2;
  while(max(boundaryDepth) < maxDepth){

    layerThickness[i] <-  layerThickness[i-1] + (layerThickness[i-1]*boundaryFactor);
    boundaryDepth[i]  <-  boundaryDepth[i-1] + layerThickness[i];
    midDepth[i]  <-  mean(c(boundaryDepth[i], boundaryDepth[i-1]));

    i = i + 1;
  }

  if(tail(boundaryDepth,maxDepth) > 1 & boundaryFactor > 1e-5){
    midDepth[length(midDepth)] <- mean(c(tail(midDepth,2)[1], maxDepth));
    layerThickness[length(layerThickness)] <- (maxDepth) - tail(boundaryDepth,2)[1];
    boundaryDepth[length(boundaryDepth)]  <- maxDepth;
  }

  boundaryDepth <- c(0,boundaryDepth)

  return(list("boundaryDepth" = boundaryDepth,
              "midDepth" = midDepth,
              "layerThickness" = layerThickness))

}

# ------------------------------------------------------------------------------
# Calculate the bulk density depth profile
# ------------------------------------------------------------------------------

calculateBulkDensity <- function(par){
  
  # The number of bulk density values provided
  r <- length(par@denityMin$bulkDensity)
  
  # If only 1 bulk density value is provided as input, a constant BD profile is created
  if(r == 1){
    bulkDensityProfile_highRes <- rep(1,par@layerNumber)
    bulkDensityProfile_highRes <- bulkDensityProfile_highRes*par@denityMin$bulkDensity
  } else { # If a depth profile of BD is provided
    lowerDepth <- par@denityMin$lowerDepth
    upperDepth <- par@denityMin$upperDepth
    midDepth <- rowMeans(cbind(upperDepth, lowerDepth))
    
    # First, the depth profile is calculated in layers of 1 cm
    bulkDensityProfile_highRes <- rep(NA,ceiling(par@totalDepth)*100)
    
    # The BD depth profile is interpolated using a spline interpolation
    BD_spline <- spline(x = midDepth, y = par@denityMin$bulkDensity, xout = seq(1,max(upperDepth)))
    splineDepth <- BD_spline$x
    splineBD <- BD_spline$y
    
    # The BD of the soil layers below the maximum depth is kept constant
    bulkDensityProfile_highRes[1:length(splineBD)] <- splineBD
    bulkDensityProfile_highRes[(length(splineBD)+1):length(bulkDensityProfile_highRes)] <- tail(splineBD,1)
    
    # The mean BD per depth layer is calculated
    bulkDensityProfile <- rep(NA,par@layerNumber)
    for(i in 1:par@layerNumber){
      upperRow <- round(par@boundaryDepth[i]*100)+1
      lowerRow <- round(par@boundaryDepth[i+1]*100)
      bulkDensityProfile[i] <- mean(bulkDensityProfile_highRes[upperRow:lowerRow])
    }
    
    # To check the results:
    # df <- data.frame("x" = bulkDensityProfile, "y" =paramValues$midDepth_layer)
    # 
    # ggplot(df, aes(x=x,y=y)) + 
    #          geom_line()
    
    return(bulkDensityProfile)
    
  } # End if
} # End function

# ------------------------------------------------------------------------------
# The Km values are converted to C%
# ------------------------------------------------------------------------------

kmDepthProf <- function(par){
 
  # The mass of soil in every layer
  par@soilMassArray = (par@bulkDensityProfile * 1000) * par@layerThickness
  
  # Rhizosphere
  par@KmB_root_R_profile_stage1 <- rep(par@KmB_root_R_stage1,par@layerNumber)
  par@KmU_Bioav_profile_stage1 <- rep(par@KmU_Bioav_stage1,par@layerNumber)

  if(par@nVeg == 2){
    par@KmB_root_R_profile_stage2 <- rep(par@KmB_root_R_stage2,par@layerNumber)
    par@KmU_Bioav_profile_stage2 <- rep(par@KmU_Bioav_stage2,par@layerNumber)
  }
  
  # Mineral soil
  par@Km_depol_M_profile_stage1 <- par@soilMassArray * par@Km_depol_M_stage1

  if(par@nVeg == 2){
    par@Km_depol_M_profile_stage2 <- par@soilMassArray * par@Km_depol_M_stage2
  }
  
  par@Km_ads_profile = par@soilMassArray * par@Km_ads

  return(par)
}

# ------------------------------------------------------------------------------
# Depth profiles of CO2 diffusivity
# ------------------------------------------------------------------------------

calculate_CO2diff <- function(par){
  
  # The depth profile of soil temperatures is calculated
  soilTempDepthProfile <- rep(1,par@layerNumber)*par@Tair

  # The depth profile of air-filled porosity is calculated, exponentially declining with depth
  totalPorosity <- 1- (par@bulkDensityProfile/2.65)
  airFilledPorosityDepthProfile <- par@airFilledPorosity_surface*exp(-(1/1)*par@midDepth_layer)
  
  # The CO2 diffusivity depth profile is calculated, according to Moldrup et al. (2000)
  stdDiffusivity <- 1.381e-5
  D0 <- (stdDiffusivity)*(101.325/par@atmPressure)*((soilTempDepthProfile+273.15)/273.15)^1.81         # Depth profile
  D0_atm <- (stdDiffusivity)*(101.325/par@atmPressure)*((par@Tair+273.15)/273.15)^1.81                     # For the atmosphere
  
  D0_atm <- D0_atm*60*60*24 # m2 per day
  
  D0_atm_year<- D0_atm*365
  
  # The effective diffusivity, corrected for air-filled pore space and moisture content (m2 per second)
  De_12C <- D0*0.66*airFilledPorosityDepthProfile*(airFilledPorosityDepthProfile/totalPorosity)^3
  De_13C <- De_12C/1.0044
  De_14C <- De_12C/1.0088
  
  # These are converted to m2 per day
  De_12C_day<- De_12C*60*60*24
  De_13C_day <- De_13C*60*60*24
  De_14C_day <- De_14C*60*60*24
  
  # And per year
  De_12C_year <- De_12C_day*365
  De_13C_year <- De_13C_day*365
  De_14C_year <- De_14C_day*365
  
  par@D0_atm_year <- D0_atm_year
  par@De_12C_year <- De_12C_year
  par@De_13C_year <- De_13C_year
  par@De_14C_year <- De_14C_year
  par@airFilledPorosityDepthProfile <- airFilledPorosityDepthProfile
  
  return(par)
}

# ------------------------------------------------------------------------------
# Creating the depth profiles of root carbon inputs
# ------------------------------------------------------------------------------

rootInputs <- function(par, vegData_stage1, vegData_stage2){
  
  # ------------------
  # Vegetation stage 1
  # ------------------
  
  # If the root profile is calculated exponentially
  if (vegData_stage1@calculateRootProfileExponentially == 1){
    
    # Write codes if needed
    
  } else if(vegData_stage1@calculateRootProfileAssymptotically){
   
    tmpProfile = 1-vegData_stage1@betaRoots^((tail(par@boundaryDepth,-1))*100)
    relProfile <- rep(NA,par@layerNumber)
    relProfile[1] <- tmpProfile[1]
    
    for(i in 2:par@layerNumber){
      relProfile[i] <- tmpProfile[i] - tmpProfile[i-1]
    } # End for
    
    # The fraction of C inputs down to the simulated depth is calculated
    fractOfTotalInputs <- tail(cumsum(relProfile),1)
    
    # The depth profile of carbon inputs is calculated
    vegData_stage1@i_bgveg_stage1 <- relProfile * vegData_stage1@i_bgveg_tot# * fractOfTotalInput
    
  } # End if
  
  # ------------------
  # Vegetation stage 2
  # ------------------
  
  if(par@nVeg == 2){

    # If the root profile is calculated exponentially
    if (vegData_stage2@calculateRootProfileExponentially == 1){

      # Write codes when needed

    } else if(vegData_stage2@calculateRootProfileAssymptotically){

      tmpProfile = 1-vegData_stage2@betaRoots^((tail(par@boundaryDepth,-1))*100)
      relProfile <- rep(NA,par@layerNumber)
      relProfile[1] <- tmpProfile[1]

      for(i in 2:par@layerNumber){
        relProfile[i] <- tmpProfile[i] - tmpProfile[i-1]
      } # End for
      
      # The fraction of C inputs down to the simulated depth is calculated
      fractOfTotalInputs <- tail(cumsum(relProfile),1)
      
      vegData_stage1@i_bgveg_stage1 <- relProfile * vegData_stage2@i_bgveg_tot# * fractOfTotalInputs

    } # End if(vegData_stage2$calculateRootProfileExponentially == 1)

  } # End if(paramValues$nVeg == 2)
  
  return(list("vegData_stage1" = vegData_stage1, "vegData_stage2" = vegData_stage2))
  
} # End function

# ------------------------------------------------------------------------------
# The annual root CO2 production depth profile
# ------------------------------------------------------------------------------

annualrootCO2Production <- function(par, vegData_stage1, vegData_stage2){
  
  # ------------------
  # Vegetation stage 1
  # ------------------
  
  relRootBiomass <- vegData_stage1@i_bgveg_stage1/sum(vegData_stage1@i_bgveg_stage1)
  annualRootCO2_stage1 <- relRootBiomass*vegData_stage1@i_bgveg_tot/2
  
  vegData_stage1@annualRootCO2_stage1 <- annualRootCO2_stage1
  
  if(par@nVeg == 2){
  
    # ------------------
    # Vegetation stage 2
    # ------------------
    
    relRootBiomass <- vegData_stage2@i_bgveg_stage2/sum(vegData_stage2@i_bgveg_stage2)
    annualRootCO2_stage2 <- relRootBiomass*vegData_stage2@i_bgveg_tot/2
    
    vegData_stage2@annualRootCO2_stage2 <- annualRootCO2_stage2
  
  } # End if(paramValues$nVeg == 2)
  
  return(list("vegData_stage1" = vegData_stage1, "vegData_stage2" = vegData_stage2))
}

# ------------------------------------------------------------------------------
# The biodiffusion coefficient depth profile
# ------------------------------------------------------------------------------

bioDiff <- function(par){
  
  # ------------------
  # Vegetation stage 1
  # ------------------
  
  DbArray_stage1 <- rep(NA,par@layerNumber)
  for(i in 1:par@layerNumber){
    DbArray_stage1[i] <- par@Db0_stage1*exp(-(1/par@Db_eFold_depth_stage1)*par@midDepth_layer[i])
  }
  
  par@DbArray_stage1 <- DbArray_stage1
  
  if(par@nVeg == 2){
    
    # ------------------
    # Vegetation stage 2
    # ------------------
    
    DbArray_stage2 <- rep(NA,par@layerNumber)
    for(i in 1:par@layerNumber){
      DbArray_stage2[i] <- par@Db0_stage2*exp(-(1/par@Db_eFold_depth_stage2)*par@midDepth_layer[i])
    }
    
    par@DbArray_stage2 <- DbArray_stage2
    
  } # End if(paramValues$nVeg == 2)
 
  return(par) 
  
} # End function

# ------------------------------------------------------------------------------
# The fraction of soil influenced by the rhizosphere
# ------------------------------------------------------------------------------
# Calculted following Finzi et al. (2018)

rhizosphereVolume <- function(par, vegData_stage1, vegData_stage2){
 
  depthInterval_tmp <- 0.01

  totalDepth_tmp <- 1 # First the distribution down to 1 m depth is calculated (as the parameters are meant to do this)
                      # At the end of the function only the volume down to the simulated depth is retained
  
  # The length of fine roots per depth layer is calculated
  depthProfile_cm <- seq(depthInterval_tmp,totalDepth_tmp,depthInterval_tmp)*100 #Unit is cm
  fineRoot_cumulDepth <- 1 - vegData_stage1@fineRootBeta^depthProfile_cm
  
  fineRoot_fractionPerDepth <- rep(NA,length(depthProfile_cm))
  fineRoot_fractionPerDepth[1] <- fineRoot_cumulDepth[1]
  fineRoot_fractionPerDepth[2:length(fineRoot_fractionPerDepth)] <- 
    fineRoot_cumulDepth[2:length(fineRoot_cumulDepth)] - fineRoot_cumulDepth[1:(length(fineRoot_cumulDepth)-1)] # Unit is a fraction (-)
     
  fineRoot_totalLengthPerDepth <- fineRoot_fractionPerDepth * vegData_stage1@fineRootLength # Unit is km/m2 per depth layer
  
  # The cumulative root diameter function is calculated
  rootDiameters <- seq(0.02,2,0.02) # Unit is mm
  # The cumulative root diameter distribution is calculated
  # First some parameters are defined (from Finzi et al., 2018)
  alpha <- 75
  gamma <- 11
  CRL = 1/(1+(alpha*exp(-gamma*rootDiameters)))
  # Proportion of CRL (cumulative root length) in each depth layer
  CRL_fractionPerDiameter <- rep(NA,length(rootDiameters))
  CRL_fractionPerDiameter[1] <- CRL[1]
  CRL_fractionPerDiameter[2:length(CRL_fractionPerDiameter)] <- 
    CRL[2:length(CRL)] - CRL[1:(length(CRL)-1)]
  
  # The exudation distance for every root diameter is calculated (fine roots exudates reach larger distances from the roots) 
  fineRoot_exudationDist <- vegData_stage1@exudationDist*exp(-vegData_stage1@exudationK*rootDiameters)
  fineRootVolume_Rhizosphere_1m <- sum(100*CRL_fractionPerDiameter*pi*((rootDiameters + fineRoot_exudationDist)/10)^2) # Unit is cm3 per m2
  
  # The volume of the rhizosphere per depth layer is calculated
  rhizosphereVolumePerLayer <- (fineRoot_totalLengthPerDepth*1000)*fineRootVolume_Rhizosphere_1m # Unit is cm3 per layer
  
  # The fraction of volume occupied by the rhizosphere in every depth layer is calculated
  volumeInCm3 <- totalDepth_tmp * (depthInterval_tmp*100)*100^2 # Down to 1 m
  
  rhizosphereVolumeFraction_vegStage1_highRes <- rhizosphereVolumePerLayer/volumeInCm3; # /10000 because this is the volume 
  
  # The average rhizosphere volume for the simulated depth layers is calculated
  rhizosphereVolumeFraction_vegStage1 <- rep(NA,par@layerNumber)
  for(i in 1:par@layerNumber){
    upperRow<- round(par@boundaryDepth[i]*100)+1
    lowerRow<- round(par@boundaryDepth[i+1]*100)
    rhizosphereVolumeFraction_vegStage1[i] = mean(rhizosphereVolumeFraction_vegStage1_highRes[upperRow:lowerRow])
  }
  
  vegData_stage1@rhizosphereVolumeFraction_vegStage1 <- rhizosphereVolumeFraction_vegStage1
  
  # df <- data.frame("x" = rhizosphereVolumeFraction_vegStage1_highRes, "y" = depthProfile_cm)
  # ggplot(df, aes(x=x, y=-y)) +
  #   geom_line()
 
  if(par@nVeg == 2){
    
    # ------------------
    # Vegetation stage 2
    # ------------------
    
    depthInterval_tmp <- 0.01
    # totalDepth_tmp <- ceiling(totalDepth*100)/100
    totalDepth_tmp <- 1 # First the distribution down to 1 m depth is calculated (as the parameters are meant to do this)
                        # At the end of the function only the volume down to the simulated depth is retained
    
    # The length of fine roots per depth layer is calculated
    depthProfile_cm <- seq(depthInterval_tmp,totalDepth_tmp,depthInterval_tmp)*100 #Unit is cm
    fineRoot_cumulDepth <- 1 - vegData_stage2$fineRootBeta^depthProfile_cm
    
    fineRoot_fractionPerDepth <- rep(NA,length(depthProfile_cm))
    fineRoot_fractionPerDepth[1] <- fineRoot_cumulDepth[1]
    fineRoot_fractionPerDepth[2:length(fineRoot_fractionPerDepth)] <- 
      fineRoot_cumulDepth[2:length(fineRoot_cumulDepth)] - fineRoot_cumulDepth[1:(length(fineRoot_cumulDepth)-1)] # Unit is a fraction (-)
    
    fineRoot_totalLengthPerDepth <- fineRoot_fractionPerDepth*vegData_stage2$fineRootLength # Unit is km/m2 per depth layer
    
    # The cumulative root diameter function is calculated
    rootDiameters <- seq(0.02,2,0.02)
    # The cumulative root diameter distribution is calculated
    # First some parameters are defined (from Finzi et al., 2018)
    alpha <- 75
    gamma <- 11
    CRL = 1/(1+(alpha*exp(-gamma*rootDiameters)));
    # Proportion of CRL (cumulative root length) in each depth layer
    CRL_fractionPerDiameter <- rep(NA,length(rootDiameters))
    CRL_fractionPerDiameter[1] <- CRL[1]
    CRL_fractionPerDiameter[2:length(CRL_fractionPerDiameter)] <- 
      CRL[2:length(CRL)] - CRL[1:(length(CRL)-1)]
    
    # The exudation distance for every root diameter is calculated (fine roots exudates reach larger distances from the roots) 
    fineRoot_exudationDist <- vegData_stage2$exudationDist*exp(-vegData_stage2$exudationK*rootDiameters)
    fineRootVolume_Rhizosphere_1m <- sum(100*CRL_fractionPerDiameter*pi*((rootDiameters + fineRoot_exudationDist)/10)^2) # Unit is cm3 per m2
    
    # The volume of the rhizosphere per depth layer is calculated
    rhizosphereVolumePerLayer <- (fineRoot_totalLengthPerDepth*1000)*fineRootVolume_Rhizosphere_1m # Unit is cm3 per layer
    
    # The fraction of volume occupied by the rhizosphere in every depth layer is calculated
    # volumeInCm3 <- totalDepth*(depthInterval_tmp*100)*100^2
    volumeInCm3 <- totalDepth_tmp * (depthInterval_tmp*100)*100^2 # Down to 1 m
    rhizosphereVolumeFraction_vegStage2_highRes <- rhizosphereVolumePerLayer/volumeInCm3; # /10000 because this is the volume 
    
    # The average rhizosphere volume for the simulated depth layers is calculated
    rhizosphereVolumeFraction_vegStage2 <- rep(NA,par@layerNumber)
    for(i in 1:par@layerNumber){
      upperRow<- round(boundaryDepth[i]*100)+1
      lowerRow<- round(boundaryDepth[i+1]*100)
      rhizosphereVolumeFraction_vegStage2[i] = mean(rhizosphereVolumeFraction_vegStage2_highRes[upperRow:lowerRow])
    }
    
    vegData_stage2@rhizosphereVolumeFraction_vegStage2 <- rhizosphereVolumeFraction_vegStage2
    
    # df <- data.frame("x" = rhizosphereVolumeFraction_vegStage1_highRes, "y" = depthProfile_cm)
    # ggplot(df, aes(x=x, y=-y)) +
    #   geom_line()
    
  } # End (par@nVeg == 2)
  
  return(list("vegData_stage1" = vegData_stage1, "vegData_stage2" = vegData_stage2))
  
} # End function

# ------------------------------------------------------------------------------
# Adsorption and desorption parameters are calculated
# ------------------------------------------------------------------------------

adsParam <- function(par, vegData_stage1, vegData_stage2){
  
  # The depth profile of CminMax is calculated, by mulitplying a C% with the soil mass per layer
  par@CminMaxProfile <- par@CminMax * par@bulkDensityProfile * 1000 * par@layerThickness
  
  # Decline of Kdes with depth
  # Following the relative volume of total soil occupied by the rhizosphere
  par@kDes_stage1 <- par@kDes_init * vegData_stage1@rhizosphereVolumeFraction_vegStage1
  
  if(par@nVeg == 2){
    par@kDes_stage2 <- par@kDes_init * vegData_stage2@rhizosphereVolumeFraction_vegStage2
  }
  
  # Decline of Kads with depth
  par@Vmax_ads_stage1 <- par@Vmax_ads# * vegData_stage1@rhizosphereVolumeFraction_vegStage1
  
  if(par@nVeg == 2){
    par@Vmax_ads_stage2 <- par@Vmax_ads * vegData_stage2@rhizosphereVolumeFraction_vegStage2
  }
  
  return(par)
}

# ------------------------------------------------------------------------------
# The depth profiles of 12CO2, 13CO2 and 14CO2 are initialized
# ------------------------------------------------------------------------------

# By running the diffusion module with root inputs as the only CO2 inputs
# (thereby neglecting inputs from mineralized SOC)

initialCO2profile <- function(par, vegData_stage1){
  
  # The initial d13C value of soil CO2 equals that of belowground vegetation
  d13C_soilCO2_init <- rep(vegData_stage1@BGveg_d13C[1], par@layerNumber)
  
  # The initial concentration of soil CO2 is 0
  soilCO2_init_12C_MolPerCubicMeter <- rep(0,par@layerNumber)
  soilCO2_init_13C_MolPerCubicMeter <- rep(0,par@layerNumber)
  soilCO2_init_14C_MolPerCubicMeter <- rep(0,par@layerNumber)
  
  # The annual CO2 production profile is converted from gram C to mole CO2 per cubic meter
  annualRootCO2Production_gramCO2 = ((vegData_stage1@annualRootCO2_stage1*1e3)/12)*44 # gram CO2 per layer
  annualRootCO2Production_Mole <- annualRootCO2Production_gramCO2/44 # Mole CO2 per soil layer
  annualRootCO2Production_MolePerCubicMeter <- annualRootCO2Production_Mole*(1/par@layerThickness)*(1/par@airFilledPorosityDepthProfile) # Mole per cubic meter
  
  # This is divided between 12CO2, 13CO2 and 14CO2
  annualRootCO2Production_12C <- annualRootCO2Production_MolePerCubicMeter /(1 + par@init_rhizo_13C_12C + par@init_rhizo_14C_12C)
  annualRootCO2Production_13C <- (annualRootCO2Production_MolePerCubicMeter / annualRootCO2Production_12C - par@init_rhizo_14C_12C - 1) * annualRootCO2Production_12C
  annualRootCO2Production_14C <- annualRootCO2Production_MolePerCubicMeter - annualRootCO2Production_12C - annualRootCO2Production_13C
  
  # The atmospheric CO2 concentration is converted to mol per cubic meter
  atmCO2Conc_MolPerCubicMeter = (par@CO2conc_atmosphere[1]*42.2963)/1e6
  
  # The concentration of 12CO2, 13CO2 and 14CO2 in the atmosphere
  atmCO2_12C_MolPerCubicMeter <- atmCO2Conc_MolPerCubicMeter / (1 + par@atm_ratio_13C_12C[1] +par@atm_ratio_14C_12C[1])
  atmCO2_13C_MolPerCubicMeter <- (atmCO2Conc_MolPerCubicMeter / atmCO2_12C_MolPerCubicMeter - par@atm_ratio_14C_12C[1] - 1) * atmCO2_12C_MolPerCubicMeter
  atmCO2_14C_MolPerCubicMeter <- atmCO2Conc_MolPerCubicMeter - atmCO2_12C_MolPerCubicMeter - atmCO2_13C_MolPerCubicMeter
  
  pars <- list("De_12C_year" = par@De_12C_year,
               "De_13C_year" = par@De_13C_year,
               "De_14C_year" = par@De_14C_year,
               "atmCO2_12C_MolPerCubicMeter" = atmCO2_12C_MolPerCubicMeter,
               "atmCO2_13C_MolPerCubicMeter" = atmCO2_13C_MolPerCubicMeter,
               "atmCO2_14C_MolPerCubicMeter" = atmCO2_14C_MolPerCubicMeter,
               "D0_atm_year_12C" = par@D0_atm_year,
               "D0_atm_year_13C" = par@D0_atm_year/1.0044,
               "D0_atm_year_14C" = par@D0_atm_year/1.0088,
               "annualRootCO2Production_12C" = annualRootCO2Production_12C,
               "annualRootCO2Production_13C" = annualRootCO2Production_13C,
               "annualRootCO2Production_14C" = annualRootCO2Production_14C
                )
  
  # The initial CO2 depth profile is calculated using and ode solver
  yini <- c("CO2_12" = soilCO2_init_12C_MolPerCubicMeter,
            "CO2_13" = soilCO2_init_13C_MolPerCubicMeter,
            "CO2_14" = soilCO2_init_14C_MolPerCubicMeter
            )
  
  times <- seq(1,1000)
  
  out <- ode(func = model_initCO2DepthProfile, y = yini, parms = pars, times = times, paramValues = par)
 
  # The output is formatted
  CO2_12 <- out[nrow(out),startsWith(colnames(out), "CO2_12")]
  CO2_13 <- out[nrow(out),startsWith(colnames(out), "CO2_13")]
  CO2_14 <- out[nrow(out),startsWith(colnames(out), "CO2_14")]
  
  # ggplot(df, aes(y = depth)) +
  #   geom_line(aes(x = d13CO2)) +
  #   geom_point(aes(x = d13CO2)) +
  #   scale_y_continuous(trans = "reverse")
  
  # df <- list("CO2_12_init" = unname(CO2_12),
  #                 "CO2_13_init" = unname(CO2_13),
  #                 "CO2_14_init" = unname(CO2_14))
  
  par@CO2_12_init <- unname(CO2_12)
  par@CO2_13_init <- unname(CO2_13)
  par@CO2_14_init <- unname(CO2_14)
  
  return(par)
  
}

# ------------------------------------------------------------------------------
# Calculate the d13C value
# ------------------------------------------------------------------------------

calculateD13C <- function(C12, C13){

  d13C <- (((C13 / C12) / 0.0112372) - 1) * 1000
    
  return(d13C)
}

# ------------------------------------------------------------------------------
# Calculate the d14C value
# ------------------------------------------------------------------------------

calculateD14C <- function(C12, C14, d13C, year){
  
  R_25 = (C14 / C12) * (0.9750 / (1 + (d13C / 1000)))^2;
  R_std = 0.95 * (1.18 * 10^(-12))*exp((year - 1950) / 8267);
  
  d14C = ((R_25 / R_std) - 1)*1000;
  
  return(d14C)
}

# ------------------------------------------------------------------------------
# The state variables are initialized
# ------------------------------------------------------------------------------

initialization <- function(par, vegData_stage1, isotopeRatios, isotopicInputs){
  
  initialStateVar <- c(
              # Litter
              "POC_12L" = NA,
              "DOC_12L" = NA,
              "MIC_12L" = NA,
              "POC_13L" = NA,
              "DOC_13L" = NA,
              "MIC_13L" = NA,
              "POC_14L" = NA,
              "DOC_14L" = NA,
              "MIC_14L" = NA,
              
              # Rhizosphere
              "POC_12R" = rep(NA, par@layerNumber),
              "BIOAV_12R" = rep(NA, par@layerNumber),
              "MIC_12R" = rep(NA, par@layerNumber),
              "CO2_12R" = rep(NA, par@layerNumber),
              "POC_13R" = rep(NA, par@layerNumber),
              "BIOAV_13R" = rep(NA, par@layerNumber),
              "MIC_13R" = rep(NA, par@layerNumber),
              "CO2_13R" = rep(NA, par@layerNumber),
              "POC_14R" = rep(NA, par@layerNumber),
              "BIOAV_14R" = rep(NA, par@layerNumber),
              "MIC_14R" = rep(NA, par@layerNumber),
              "CO2_14R" = rep(NA, par@layerNumber),
              
              # Mineral soil
              "POLY_12B" = rep(NA, par@layerNumber),
              "MIN_12B" = rep(NA, par@layerNumber),
              "MIC_12B" = rep(NA, par@layerNumber),
              "CO2_12B" = rep(NA, par@layerNumber),
              "SURF_12B" = rep(NA, par@layerNumber),
              "POLY_13B" = rep(NA, par@layerNumber),
              "MIN_13B" = rep(NA, par@layerNumber),
              "MIC_13B" = rep(NA, par@layerNumber),
              "CO2_13B" = rep(NA, par@layerNumber),
              "SURF_13B" = rep(NA, par@layerNumber),
              "POLY_14B" = rep(NA, par@layerNumber),
              "MIN_14B" = rep(NA, par@layerNumber),
              "MIC_14B" = rep(NA, par@layerNumber),
              "CO2_14B" = rep(NA, par@layerNumber),
              "SURF_14B" = rep(NA, par@layerNumber))
  
  # A small amount of carbon is added
  
  # Microbes in the litter layer
  Cmic_litter_init <- (vegData_stage1@i_agv/20)
  initialStateVar["MIC_13L"] <- Cmic_litter_init / (1 + ( 1 / isotopeRatios$BGveg_ratio_13C_12C_deadRoots[1]))
  initialStateVar["MIC_12L"] <- initialStateVar["MIC_13L"] / isotopeRatios$BGveg_ratio_13C_12C_deadRoots[1]
  initialStateVar["MIC_14L"] <- Cmic_litter_init / (1 + (1 / isotopeRatios$BGveg_ratio_14C_12C_deadRoots[1]))
  
  # Microbes in the rhizosphere
  Cmic_init <- (par@soilMassArray * 1e-6 * vegData_stage1@rhizosphereVolumeFraction_vegStage1)
  initialStateVar[startsWith(names(initialStateVar), "MIC_13R")] <- Cmic_init / (1 + (1 / isotopeRatios$BGveg_ratio_13C_12C_deadRoots[1]))
  initialStateVar[startsWith(names(initialStateVar), "MIC_12R")] <- initialStateVar[startsWith(names(initialStateVar), "MIC_13R")] / isotopeRatios$BGveg_ratio_13C_12C_deadRoots[1]
  initialStateVar[startsWith(names(initialStateVar), "MIC_14R")] <- Cmic_init / (1 + (1 / isotopeRatios$BGveg_ratio_14C_12C_deadRoots[1]))
  
  # Microbes in the bulk soil
  initialStateVar[startsWith(names(initialStateVar), "MIC_13B")] <- Cmic_init / (1 + (1 / isotopeRatios$BGveg_ratio_13C_12C_deadRoots[1]))
  initialStateVar[startsWith(names(initialStateVar), "MIC_12B")] <- initialStateVar[startsWith(names(initialStateVar), "MIC_13B")] / isotopeRatios$BGveg_ratio_13C_12C_deadRoots[1]
  initialStateVar[startsWith(names(initialStateVar), "MIC_14B")] <- Cmic_init / (1 + (1 / isotopeRatios$BGveg_ratio_14C_12C_deadRoots[1]))
  
  # POC in the litter layer
  initialStateVar["POC_12L"] <- isotopicInputs$AGveg_12C_inputs[1] * (1 - par@L_fractionLeachable)
  initialStateVar["POC_13L"] <- isotopicInputs$AGveg_13C_inputs[1] * (1 - par@L_fractionLeachable)
  initialStateVar["POC_14L"] <- isotopicInputs$AGveg_14C_inputs[1] * (1 - par@L_fractionLeachable)
  
  # DOC in the litter layer
  initialStateVar["DOC_12L"] <- isotopicInputs$AGveg_12C_inputs[1] * par@L_fractionLeachable
  initialStateVar["DOC_13L"] <- isotopicInputs$AGveg_13C_inputs[1] * par@L_fractionLeachable
  initialStateVar["DOC_14L"] <- isotopicInputs$AGveg_14C_inputs[1] * par@L_fractionLeachable
  
  # Mineral carbon in the bulk soil
  C_init = ((par@soilMassArray * vegData_stage1@rhizosphereVolumeFraction_vegStage1) / 1000)
  initialStateVar[startsWith(names(initialStateVar), "MIN_13B")] <- C_init / (1 + (1 / isotopeRatios$BGveg_ratio_13C_12C_deadRoots[1]))
  initialStateVar[startsWith(names(initialStateVar), "MIN_12B")] <- initialStateVar[startsWith(names(initialStateVar), "MIN_13B")] / isotopeRatios$BGveg_ratio_13C_12C_deadRoots[1]
  initialStateVar[startsWith(names(initialStateVar), "MIN_14B")] <- C_init / (1 + (1 / isotopeRatios$BGveg_ratio_14C_12C_deadRoots[1]))
  
  # Polymeric carbon in the bulk soil
  initialStateVar[startsWith(names(initialStateVar), "POLY_13B")] <- initialStateVar[startsWith(names(initialStateVar), "MIN_13B")]
  initialStateVar[startsWith(names(initialStateVar), "POLY_12B")] <- initialStateVar[startsWith(names(initialStateVar), "MIN_12B")]
  initialStateVar[startsWith(names(initialStateVar), "POLY_14B")] <- initialStateVar[startsWith(names(initialStateVar), "MIN_14B")]
  
  # POC in the rhizosphere
  initialStateVar[startsWith(names(initialStateVar), "POC_12R")] <- isotopicInputs$BGveg_12C_inputs_deadRoots[,1]
  initialStateVar[startsWith(names(initialStateVar), "POC_13R")] <- isotopicInputs$BGveg_13C_inputs_deadRoots[,1]
  initialStateVar[startsWith(names(initialStateVar), "POC_14R")] <- isotopicInputs$BGveg_14C_inputs_deadRoots[,1]
  
  # BIOAV in the rhizosphere
  initialStateVar[startsWith(names(initialStateVar), "BIOAV_12R")] <- isotopicInputs$BGveg_12C_inputs_rhizodeposition[,1]
  initialStateVar[startsWith(names(initialStateVar), "BIOAV_13R")] <- isotopicInputs$BGveg_13C_inputs_rhizodeposition[,1]
  initialStateVar[startsWith(names(initialStateVar), "BIOAV_14R")] <- isotopicInputs$BGveg_14C_inputs_rhizodeposition[,1]
  
  # CO2 in the rhizosphere
  initialStateVar[startsWith(names(initialStateVar), "CO2_12R")] <- par@CO2_12_init
  initialStateVar[startsWith(names(initialStateVar), "CO2_13R")] <- par@CO2_13_init
  initialStateVar[startsWith(names(initialStateVar), "CO2_14R")] <- par@CO2_14_init
  
  # CO2 in the bulk soil
  initialStateVar[startsWith(names(initialStateVar), "CO2_12B")] <- par@CO2_12_init#rep(0,par@layerNumber)
  initialStateVar[startsWith(names(initialStateVar), "CO2_13B")] <- par@CO2_13_init#rep(0,par@layerNumber)
  initialStateVar[startsWith(names(initialStateVar), "CO2_14B")] <- par@CO2_14_init#rep(0,par@layerNumber)
  
  # The maximum amount of mineral-associated C is added
  initialStateVar[startsWith(names(initialStateVar), "SURF_13B")] <- par@CminMaxProfile / (1 + (1 / isotopeRatios$BGveg_ratio_13C_12C_deadRoots[1]))
  initialStateVar[startsWith(names(initialStateVar), "SURF_12B")] <- initialStateVar[startsWith(names(initialStateVar), "SURF_13B")] / isotopeRatios$BGveg_ratio_13C_12C_deadRoots[1]
  initialStateVar[startsWith(names(initialStateVar), "SURF_14B")] <- par@CminMaxProfile / (1 + (1 / isotopeRatios$BGveg_ratio_14C_12C_deadRoots[1]))

  return(initialStateVar)
  
}

# ------------------------------------------------------------------------------
# A function that assigns the correct parameters per vegetation stage
# ------------------------------------------------------------------------------

stageParameters <- function(par, vegData_stage1, vegData_stage2, stage){
  
  if(stage == 1){
    
    p <- list("Kdes" = par@kDes_stage1,
              "Vmax_ads" = par@Vmax_ads_stage1,
              "Km_ads" = par@Km_ads_profile,
              "kDes" = par@kDes_stage1,
              "Db0" = par@Db0_stage1,
              "Db_eFold_depth" = par@Db_eFold_depth_stage1,
              "bioturbFract" = par@bioturbFract_stage1,
              "advectionRate_polyC" = par@advectionRate_polyC_stage1,
              "Vmax_L" = par@Vmax_L,
              "VmaxD_POC_L" = par@VmaxD_POC_L_stage1,
              "KmB_POC_L" = par@KmB_POC_L_stage1,
              "VmaxD_doc_L" = par@VmaxD_doc_L_stage1,
              "KmB_DOC_L" = par@KmB_DOC_L_stage1,
              "K_mic_L" = par@K_mic_L,
              "K_mic_R" = par@K_mic_R,
              "K_mic_B" = par@K_mic_B,
              "bioToPoly" = par@bioToPoly,
              "leachFract" = par@leachFract_stage1,
              "VmaxD_root_R" = par@VmaxD_root_R_stage1,
              "KmB_root_R" = par@KmB_root_R_stage1,
              "VmaxU_BioAv" = par@VmaxU_BioAv_stage1,
              "KmU_Bioav_profile" = par@KmU_Bioav_profile_stage1,
              "VmaxD_M" = par@VmaxD_M_stage1,
              "Km_depol_M" = par@Km_depol_M_stage1,
              "KmB_root_R_profile" = par@KmB_root_R_profile_stage1,
              "Km_depol_M_profile" = par@Km_depol_M_profile_stage1,
              "DbArray" = par@DbArray_stage1,
              "annualRootCO2" = vegData_stage1@annualRootCO2_stage1,
              "rhizosphereVolumeFraction" = vegData_stage1@rhizosphereVolumeFraction_vegStage1,
              
              # Other parameters needed in SOILcarb
              "CUE_L" = par@CUE_L,
              "CUE_R" = par@CUE_R,
              "CUE_B" = par@CUE_B,
              "fractIncorpCO2_litter" = par@fractIncorpCO2_litter,
              "fractIncorpCO2_soil" = par@fractIncorpCO2_soil,
              "solubleMicFraction" = par@solubleMicFraction,
              "beta" = par@beta,
              "L_fractionLeachable" = par@L_fractionLeachable,
              "Fdoc_rhizo" = par@Fdoc_rhizo,
              "includeBioturbation" = par@includeBioturbation,
              "includeAdvection" = par@includeAdvection,
              "d13C_atmosphere" = par@d13C_atmosphere,
              "d14C_atmosphere" = par@d14C_atmosphere,
              "layerNumber" = par@layerNumber,
              "CminMaxProfile" = par@CminMaxProfile,
              "airFilledPorosityDepthProfile" = par@airFilledPorosityDepthProfile,
              "layerThickness" = par@layerThickness,
              "midDepth_layer" = par@midDepth_layer,
              "boundaryDepth" = par@boundaryDepth
    )
    
    
  }else if(stage == 2){
    
    p <- list("Kdes" = par@kDes_stage2,
              "Vmax_ads" = par@Vmax_ads_stage2,
              "Km_ads" = par@Km_ads_profile,
              "Db0" = par@Db0_stage2,
              "Db_eFold_depth" = par@Db_eFold_depth_stage2,
              "bioturbFract" = par@bioturbFract_stage2,
              "advectionRate_polyC" = par@advectionRate_polyC_stage2,
              "VmaxD_POC_L" = par@VmaxD_POC_L_stage2,
              "KmB_POC_L" = par@KmB_POC_L_stage2,
              "VmaxD_doc_L" = par@VmaxD_doc_L_stage2,
              "KmB_doc_L" = par@KmB_doc_L_stage2,
              "bioToPoly" = par@bioToPoly,
              "leachFract" = par@leachFract_stage2,
              "VmaxD_root_R" = par@VmaxD_root_R_stage2,
              "KmB_root_R" = par@KmB_root_R_stage2,
              "VmaxU_BioAv" = par@VmaxU_BioAv_stage2,
              "KmU_Bioav_profile" = par@KmU_Bioav_profile_stage2,
              "VmaxD_M" = par@VmaxD_M_stage2,
              "Km_depol_M" = par@Km_depol_M_stage2,
              "KmB_root_R_profile" = par@KmB_root_R_profile_stage2,
              "Km_depol_M_profile" = par@Km_depol_M_profile_stage2,
              "DbArray" = par@DbArray_stage2,
              "annualRootCO2" = vegData_stage2@annualRootCO2_stage2,
              "rhizosphereVolumeFraction" = vegData_stage2@rhizosphereVolumeFraction_vegStage1,
              
              # Other parameters needed in SOILcarb
              "CUE_L" = par@CUE_L,
              "CUE_R" = par@CUE_R,
              "CUE_B" = par@CUE_B,
              "fractIncorpCO2_litter" = par@fractIncorpCO2_litter,
              "fractIncorpCO2_soil" = par@fractIncorpCO2_soil,
              "solubleMicFraction" = par@solubleMicFraction,
              "beta" = par@beta,
              "kMic_L" = par@kMic_L,
              "kMic_R" = par@kMic_R,
              "kMic_B" = par@kMic_B,
              "L_fractionLeachable" = par@L_fractionLeachable,
              "Fdoc_rhizo" = par@Fdoc_rhizo,
              "includeBioturbation" = par@includeBioturbation,
              "includeAdvection" = par@includeAdvection,
              "d13C_atmosphere" = par@d13C_atmosphere,
              "d14C_atmosphere" = par@d14C_atmosphere,
              "layerNumber" = par@layerNumber,
              "CminMaxProfile" = par@CminMaxProfile,
              "airFilledPorosityDepthProfile" = par@airFilledPorosityDepthProfile,
              "layerThickness" = par@layerThickness,
              "midDepth_layer" = par@midDepth_layer,
              "boundaryDepth" = par@boundaryDepth
    )
    
  }
  
  return(p)
  
}

# ------------------------------------------------------------------------------
# The output of the SOILcarb function is formatted
# ------------------------------------------------------------------------------

formatOutput <- function(out, year, vegData_stage1){
  # Out is a deSolve structure
  
  # --------------------------------------------------
  # The results for the different pools are structured
  # --------------------------------------------------
  
  res <- list(
    
    # ----------------
    # The litter layer
    # ----------------
    
    "litter_12C" = out[,"POC_12L"] + out[,"MIC_12L"] + out[,"DOC_12L"],
    "litter_13C" = out[,"POC_13L"] + out[,"MIC_13L"] + out[,"DOC_13L"],
    "litter_14C" = out[,"POC_14L"] + out[,"MIC_14L"] + out[,"DOC_14L"],
    
    # ---------------
    # The soil layers
    # ---------------
    
    "soil_12C" = unname(getOutput(out,"POC_12R") + getOutput(out,"BIOAV_12R") + getOutput(out,"MIC_12R") + 
                      getOutput(out,"POLY_12B") + getOutput(out,"MIC_12B") + getOutput(out,"MIN_12B")),
    "soil_13C" = unname(getOutput(out,"POC_13R") + getOutput(out,"BIOAV_13R") + getOutput(out,"MIC_13R") + 
                          getOutput(out,"POLY_13B") + getOutput(out,"MIC_13B") + getOutput(out,"MIN_13B")),
    "soil_14C" = unname(getOutput(out,"POC_14R") + getOutput(out,"BIOAV_14R") + getOutput(out,"MIC_14R") + 
                          getOutput(out,"POLY_14B") + getOutput(out,"MIC_14B") + getOutput(out,"MIN_14B")),
    
    # ----------------------------------
    # The total amount of carbon per pool
    # ----------------------------------
    
    "POC_L" = out[,"POC_12L"] + out[,"POC_13L"] + out[,"POC_14L"],
    "DOC_L" = out[,"DOC_12L"] + out[,"DOC_13L"] + out[,"DOC_14L"],
    "MIC_L" = out[,"MIC_12L"] + out[,"MIC_13L"] + out[,"MIC_14L"],
    
    "POC_R" = unname(getOutput(out,"POC_12R") + getOutput(out,"POC_13R") + getOutput(out,"POC_14R")),
    "BIOAV_R" = unname(getOutput(out,"BIOAV_12R") + getOutput(out,"BIOAV_13R") + getOutput(out,"BIOAV_14R")),
    "MIC_R" = unname(getOutput(out,"MIC_12R") + getOutput(out,"MIC_13R") + getOutput(out,"MIC_14R")),
    
    "POLY_B" = unname(getOutput(out,"POLY_12B") + getOutput(out,"POLY_13B") + getOutput(out,"POLY_14B")),
    "MIN_B" = unname(getOutput(out,"MIN_12B") + getOutput(out,"MIN_13B") + getOutput(out,"MIN_14B")),
    "MIC_B" = unname(getOutput(out,"MIC_12B") + getOutput(out,"MIC_13B") + getOutput(out,"MIC_14B")),
    
    # ---
    # CO2
    # ---
    
    "CO2_R" = unname(getOutput(out,"CO2_12R") + getOutput(out,"CO2_13R") + getOutput(out,"CO2_14R")),
    "CO2_B" = unname(getOutput(out,"CO2_12B") + getOutput(out,"CO2_13B") + getOutput(out,"CO2_14B")),
    
    "CO2_12" = unname(getOutput(out,"CO2_12R") + getOutput(out,"CO2_12B")),
    "CO2_13" = unname(getOutput(out,"CO2_13R") + getOutput(out,"CO2_13B")),
    "CO2_14" = unname(getOutput(out,"CO2_14R") + getOutput(out,"CO2_14B"))
  )
  
  # The turnover rate of microbes (for 12C) in the rhizosphere and bulk soil
  # Rhizosphere
  rowNum <- dim(res[["MIC_R"]])[1]
  Ctot_R <- res[["MIC_R"]][rowNum,] + res[["BIOAV_R"]][rowNum,] + res[["POC_R"]][rowNum,]
  K_R <- par@K_mic_R * res[["BIOAV_R"]][rowNum,]
  MicToBioAv_R <- res[["MIC_R"]][rowNum,] / res[["BIOAV_R"]][rowNum,]
  uptake_R <- par@VmaxU_BioAv_stage1 * res[["BIOAV_R"]][rowNum,] * (MicToBioAv_R / (par@KmU_Bioav_profile_stage1 + MicToBioAv_R))
  growthRate_R <- (uptake_R * par@CUE_R) / res[["MIC_R"]][rowNum,]
  death_R <- (growthRate_R * res[["MIC_R"]][rowNum,]^2) / K_R # Abs. amount of death microbes per year
  TT_rhizo_yr <- res[["MIC_R"]][rowNum,] / death_R # Tunrover time = stock / flux [year]
  TT_rhizo_day <- TT_rhizo_yr * 365 # Turnover time converted to days
  
  # Bulk soil
  Ctot_B <- res[["MIC_B"]][rowNum,] + res[["POLY_B"]][rowNum,] + res[["MIN_B"]][rowNum,]
  K_B <- par@K_mic_B * Ctot_B
  SURF_12B <- unname(getOutput(out,"SURF_12B"))[rowNum,]
  surf <- SURF_12B * vegData_stage1@rhizosphereVolumeFraction_vegStage1
  uptake_B <- (par@VmaxD_M_stage1 * res[["POLY_B"]][rowNum,] * res[["MIC_B"]][rowNum,]) / ((par@Km_depol_M_stage1 * (1 + (surf/par@Km_ads_profile) + (res[["MIC_B"]][rowNum,]/par@Km_depol_M_stage1))) + res[["POLY_B"]][rowNum,])
  growthRate_B <- (uptake_B * par@CUE_B) / res[["MIC_B"]][rowNum,]
  death_B <- (growthRate_B * res[["MIC_B"]][rowNum,]^2) / K_B # Abs. amount of death microbes per year
  TT_bulk_yr <- res[["MIC_B"]][rowNum,] / death_B # Tunrover time = stock / flux [year]
  TT_bulk_day <- TT_bulk_yr * 365 # Turnover time converted to days
  
  # Results for total carbon
  res <- c(res,
           "C_litter" = list(res$litter_12C + res$litter_13C + res$litter_14C),
           "C_soil" = list(res$soil_12C + res$soil_13C + res$soil_14C),
           "C_rhizo" = list(res$POC_R + res$BIOAV_R + res$MIC_R),
           "C_bulk" = list(res$POLY_B + res$MIN_B + res$MIC_B),
           "TT_rhizo_day" = list(TT_rhizo_day),
           "TT_bulk_day" = list(TT_bulk_day)
           )
  
  # -----------------------------------------
  # The d13C value of the pools is calculated
  # -----------------------------------------
  
  d13C <- list(
    
    # ------------
    # Litter layer
    # ------------
    
    "POC_L" = calculateD13C(out[,"POC_12L"], out[,"POC_13L"]),
    "DOC_L" = calculateD13C(out[,"DOC_12L"], out[,"DOC_13L"]),
    "MIC_L" = calculateD13C(out[,"MIC_12L"], out[,"MIC_13L"]),
    "C_L" = calculateD13C(res$litter_12C, res$litter_13C),
    
    # -----------
    # Rhizosphere
    # -----------
    
    "POC_R" = unname(calculateD13C(getOutput(out,"POC_12R"),getOutput(out,"POC_13R"))),
    "BIOAV_R" = unname(calculateD13C(getOutput(out,"BIOAV_12R"),getOutput(out,"BIOAV_13R"))),
    "MIC_R" = unname(calculateD13C(getOutput(out,"MIC_12R"),getOutput(out,"MIC_13R"))),

    # ---------
    # Bulk soil
    # ---------
    
    "POLY_B" = unname(calculateD13C(getOutput(out,"POLY_12B"),getOutput(out,"POLY_13B"))),
    "MIN_B" = unname(calculateD13C(getOutput(out,"MIN_12B"),getOutput(out,"MIN_13B"))),
    "MIC_B" = unname(calculateD13C(getOutput(out,"MIC_12B"),getOutput(out,"MIC_13B"))),
    
    # ---
    # CO2
    # ---
    
    "CO2_R" = unname(calculateD13C(getOutput(out,"CO2_12R"),getOutput(out,"CO2_13R"))),
    "CO2_B" = unname(calculateD13C(getOutput(out,"CO2_12B"),getOutput(out,"CO2_13B"))),
    "CO2_tot" = unname(calculateD13C(getOutput(out,"CO2_12R") + getOutput(out,"CO2_12B"), getOutput(out,"CO2_13R") + getOutput(out,"CO2_13B"))),
    
    # -----------
    # Soil layers
    # -----------
    "soilC" = unname(calculateD13C(res$soil_12C,res$soil_13C))
    
  )
  
  # The d13C of the bulk soil is added
  d13C <- c(d13C,
            C_bulk = list(((d13C$POLY_B * res$POLY_B) + (d13C$MIN_B * res$MIN_B) + (d13C$MIC_B * res$MIC_B)) / res$C_bulk)
              )
  
  # -----------------------------------------
  # The d14C value of the pools is calculated
  # -----------------------------------------
  
  d14C <- list(
    
    # ------------
    # Litter layer
    # ------------
    
    "POC_L" = calculateD14C(out[,"POC_12L"], out[,"POC_14L"], d13C$"POC_L", year),
    "DOC_L" = calculateD14C(out[,"DOC_12L"], out[,"DOC_14L"], d13C$"DOC_L", year),
    "MIC_L" = calculateD14C(out[,"MIC_12L"], out[,"MIC_14L"], d13C$"MIC_L", year),
    "C_L" = calculateD14C(res$litter_12C, res$litter_14C, d13C$"C_L", year),
    
    # -----------
    # Rhizosphere
    # -----------
    
    "POC_R" = unname(calculateD14C(getOutput(out,"POC_12R"),getOutput(out,"POC_14R"), d13C$"POC_R", year)),
    "BIOAV_R" = unname(calculateD14C(getOutput(out,"BIOAV_12R"),getOutput(out,"BIOAV_14R"), d13C$"BIOAV_R", year)),
    "MIC_R" = unname(calculateD14C(getOutput(out,"MIC_12R"),getOutput(out,"MIC_14R"), d13C$"MIC_R", year)),
    
    # ---------
    # Bulk soil
    # ---------
    
    "POLY_B" = unname(calculateD14C(getOutput(out,"POLY_12B"),getOutput(out,"POLY_14B"), d13C$"POLY_B", year)),
    "MIN_B" = unname(calculateD14C(getOutput(out,"MIN_12B"),getOutput(out,"MIN_14B"), d13C$"MIN_B", year)),
    "MIC_B" = unname(calculateD14C(getOutput(out,"MIC_12B"),getOutput(out,"MIC_14B"), d13C$"MIC_B", year)),
    
    # ---
    # CO2
    # ---
    
    "CO2_R" = unname(calculateD14C(getOutput(out,"CO2_12R"),getOutput(out,"CO2_14R"), d13C$"CO2_R", year)),
    "CO2_B" = unname(calculateD14C(getOutput(out,"CO2_12B"),getOutput(out,"CO2_14B"), d13C$"CO2_B", year)),
    "CO2_tot" = unname(calculateD14C(getOutput(out,"CO2_12R") + getOutput(out,"CO2_12B"), getOutput(out,"CO2_14R") + getOutput(out,"CO2_14B"), d13C$"CO2_tot", year)),
    
    # -----------
    # Soil layers
    # -----------
    
    "soilC" = unname(calculateD14C(res$soil_12C,res$soil_14C, d13C$"soilC", year))
    
  )
  
  # The d13C of the bulk soil is added
  d14C <- c(d14C,
            C_bulk = list(((d14C$POLY_B * res$POLY_B) + (d14C$MIN_B * res$MIN_B) + (d14C$MIC_B * res$MIC_B)) / res$C_bulk)
  )
  
  return(list("res" = res,
              "d13C" = d13C,
              "d14C" = d14C))
 
}


# ------------------------------------------------------------------------------
# A function to get the output from the deSolve structure
# ------------------------------------------------------------------------------

getOutput <- function(out,field){
  # Out is the deSolve output structure
  # field is the name of the variable for which the data needs to be isolated
  
  x <- out[,startsWith(colnames(out), field)]
  return(x)
  
}