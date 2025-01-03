# ------------------------------------
# A function to create the input class
# ------------------------------------

createInputClass <- function(){
  
  setClass("inputObject", slots=list(
    
    # ------------------
    # General parameters
    # ------------------
    
    CUE_L = "numeric", # Carbon use efficiency of microbes in the litter layer [-]
    CUE_R = "numeric", # Carbon use efficiency of microbes in the rhizosphere [-]
    CUE_B = "numeric", # Carbon use efficiency of microbes in the bulk soil [-]
    Db0_stage1 = "numeric", # Biodiffusion coefficient at the soil surface [m2/yr]
    Db_eFold_depth_stage1 = "numeric", # e-folding depth for the biodiffusion coefficient [m]
    denityMin = "data.frame", # Soil bulk density [g/cm3]
    advectionRate_polyC_stage1 = "numeric", # Advection velocity of DOC in the soil [m/yr]
    advectionRate_polyC_stage2 = "numeric", # Advection velocity of DOC in the soil [m/yr]
    fractIncorpCO2_litter = "numeric", # Fraction of microbial C uptake as CO2 in the litter layer [-]
    fractIncorpCO2_soil = "numeric", # Fraction of microbial C uptake as CO2 in the soil [-]
    solubleMicFraction = "numeric", # Soluble portion of microbial necromass [-]
    beta = "numeric", # Extinction coefficient for root depth profile [-]
    Db0_stage2 = "numeric", # Biodiffusion coefficient at the soil surface [m2/yr]
    Db_eFold_depth_stage2 = "numeric", # e-folding depth for the biodiffusion coefficient [m]
    CO2conc_fractionationFactor = "numeric", # Change in fractionation against 13C by plants per unit change in atmospheric CO2 concentration [per mil/ppm]
    diff_leaf_rhizoDeposits = "numeric", # The difference in δ13C values between leaves and rhizodeposits [per mil]
    
    DbArray_stage1 = "numeric", # A varibale to store the depth profile of biodiffusivity in
    
    # -------------------------------
    # Parameters for the litter layer
    # -------------------------------
    
    Vmax_L = "numeric", # Max. rate of depolimerisation of POC and DOC (same value) [yr-1]
    K_mic_L = "numeric", # Carrying capactiy as a % of litter C
    VmaxD_POC_L_stage1 = "numeric",  # Not used in the current version, Vmax_L is used instead
    KmB_POC_L_stage1 = "numeric", # Affinity for litter POC depolymerisation [kg C m-2]
    VmaxD_doc_L_stage1 = "numeric", # Not used in the current version, Vmax_L is used instead
    KmB_DOC_L_stage1 = "numeric", # Affinity for litter DOC depolymerisation [kg C m-2]
    L_fractionLeachable = "numeric", # Fraction leachable carbon of total litter inputs  [-]
    leachFract_stage1 = "numeric", # Portion of annually leached DOC from the litter layer [-]
    bioturbFract_stage1 = "numeric", # Portion of annually bioturbated POC from the litter layer [-]
    
    VmaxD_POC_L_stage2 = "numeric", # Same as above, but for the second vegetation stage (not used in this version)
    KmB_POC_L_stage2 = "numeric",
    maxD_doc_L_stage2 = "numeric",
    KmB_doc_L_stage2 = "numeric",
    leachFract_stage2 = "numeric",
    bioturbFract_stage2 = "numeric",
    
    # ------------------------------
    # Parameters for the rhizosphere
    # ------------------------------
    
    K_mic_R = "numeric", # Carrying capactiy as a % of C
    Fdoc_rhizo = "numeric", # Portion of total belowground carbon inputs as rhizodeposits [-]
    
    VmaxD_root_R_stage1 = "numeric", # Maximum rate of rhizosphere POC depolimerisation [yr-1]
    KmB_root_R_stage1 = "numeric", # Half-saturation constant for rhizosphere POC depolimerisation [ratio of MIC to POC]
    
    VmaxU_BioAv_stage1 = "numeric", # Maximum rate of carbon uptake by rhizosphere microbes [yr-1]
    KmU_Bioav_stage1 = "numeric", # Half-saturation constant for C uptake by rhizosphere microbes [MIC to BIOAV]
    
    KmB_root_R_profile_stage1 = "numeric", # To store calculated depth profiles
    KmU_Bioav_profile_stage1 = "numeric", # To store calculated depth profiles
    
    bioToPoly = "numeric", # The portion of bio-available C transferred to the soil DOC pool [yr-1]
    
    VmaxD_root_R_stage2 = "numeric", # Same as above, but for the second vegetation stage (not used in this version)
    KmB_root_R_stage2 = "numeric",
    KmB_root_R_profile_stage2 = "numeric",
    VmaxU_BioAv_stage2 = "numeric",
    KmU_Bioav_stage2 = "numeric",
    KmU_Bioav_profile_stage2 = "numeric",
    
    # ----------------------------
    # Parameters for the bulk soil
    # ----------------------------
    
    VmaxD_M_stage1 = "numeric", # Maximum rate of carbon depolimerisation and uptake in the bulk soil [yr-1]
    Km_depol_M_stage1 = "numeric", # Affinity parameter for carbon depolimerisation and uptake by microbes in the bulk soil [!!! here defined as a % of total soil mass, equivalent to %SOC]
    Vmax_ads = "numeric", # Maximum rate of adsorption of DOC onto mineral surfaces [yr-1]
    Km_ads = "numeric", # Affinity parameter for adsorption of soil DOC [!!! here defined as a % of total soil mass, equivalent to %SOC]
    kDes_init = "numeric", # Rate of desorption of mineral-associated carbon [yr-1]
    CminMax = "numeric", # Maximum amount of mineral-associated carbon [kg C m-2 per depth layer; calculated following Georgiou et al. 2022, NComm]
    Km_depol_M_profile_stage1 = "numeric", # To store the calculated depth profile of Km_depol in (after accounting for layer thickness)
    Km_ads_profile = "numeric", # To store the calculated depth profile of Km_ads in (after accounting for layer thickness)
    CminMaxProfile = "numeric", # To store the calculated depth profile of CminMax in (after accounting for layer thickness)
    kDes_stage1 = "numeric", # To store Kdes for the 1st vegetation stage
    Vmax_ads_stage1 = "numeric", # To store Vmax_ads for the 1st vegetation stage
    K_mic_B = "numeric", # Carrying capacity of microbes in the bulk soil, as a portion of C in the bulk soil [0 - 1]
    
    VmaxD_M_stage2 = "numeric", # Same as above, but for the second vegetation stage (not used in this version)
    Km_depol_M_stage2 = "numeric",
    Km_depol_M_profile_stage2 = "numeric",
    kDes_stage2 = "numeric",
    Vmax_ads_stage2 = "numeric",
    
    # -------------------------------------------------------
    # Parameters for the calculation of the CO2 depth profile
    # -------------------------------------------------------
    
    Tair = "numeric", # The average soil air temperature [k]
    airFilledPorosity_surface = "numeric", # Air-filled porosity at the soil surface [m3/m3]
    atmPressure = "numeric", # The atmospheric pressure [pascals; in the manuscript this was atm]
    dailyRootCO2Production = "numeric", # The daily amount of root CO2 production [kg CO2-C m-3 d-1]
    D0_atm_year = "numeric", # Needed in the calculations
    De_12C_year = "numeric", # Needed in the calculations
    De_13C_year = "numeric", # Needed in the calculations
    De_14C_year = "numeric", # Needed in the calculations
    airFilledPorosityDepthProfile = "numeric", # Needed in the calculations
    
    # --------------------
    # Technical parameters
    # --------------------
    
    nVeg = "numeric", # The number of simulated vegetation stages (each stage can have different vegetation characteristics)
    yearsStage1 = "numeric", # The number of years of vegetation stage 1
    yearsStage2 = "numeric", # The number of years of vegetation stage 2
    lastYearStage2 = "numeric", # The last calendar year of vegetation stage 2
    layerNumber = "numeric", # The number of soil depth layers (calculated)
    totalDepth = "numeric", # The total simulated soil depth (calculated)
    layerThickness = "numeric", # The thickness of every simulated layer (calculated)
    midDepth_layer = "numeric", # The mid depth of every simulated layer (calculated)
    lastYearStage1 = "numeric", # The last calendar year of vegetation stage 1 
    Include_Suess_effect = "numeric", # If 1, the Suess effect is simulated (if this data is provided)
    includeBioturbation = "numeric", # if 1, bioturbation is simulated
    includeAdvection = "numeric", # if 1, downward advection is simulated
    lastRealSimulationYear = "numeric", # The last calendar year of the simulations (necessary to get the right isotope inputs)
    numberOfSimulationYears = "numeric", # The number of simulated years (calculated)
    times = "numeric", # Necessary for the deSolve solver
    calibMode = "numeric", # The type of calibration mode
    evalDepth = "numeric", # The depth down to which the model has to be evaluated in calibration mode
    parallelCalib = "numeric", # If 1, the calibration is run in parallel mode
    counter_loadIsotopes = "numeric", # Necessary for the calculations
    runSeries = "numeric", # If 1, the user provided a list of parameters sets, and the model is run with all of these (in parallel)
    normalRun = "numeric", # In case the provided parameters are used for the run
    SAFE_sensitivity_isotopes = "numeric", # If 1, a SAFE sensitivity analysis is run for the isotopes
    SAFE_sensitivity_allParameters = "numeric", # If 1, a SAFE sensitivity analysis is run for all parameters
    litterCalibration = "numeric", # If 1, the parameters for only the litter layer are calibrated
    soilCalibration = "numeric", # If 1, the parameters for only the soil are calibrated
    useCalibrationResults = "numeric", # If 1, the output of the DE algorithm is used to run the model (to check the optimized model with the lowest error)
    site = "character", # The name of the site (in case data for multiple sites is provided)
    boundaryDepth = "numeric", # Calculated by the model (boundary depths of the simulated layers)
    bulkDensityProfile = "numeric", # Depth profile of the bulk density
    soilMassArray = "numeric", # Depth profile of the soil mass per simulated layer
    calibCheck = "data.frame", # To decide which calibration scenario is run
    
    # -----------------------------------------
    # Variables containing isotopic information
    # -----------------------------------------
    
    d13C_AGveg = "numeric", # The d13C value of aboveground vegetation in the last simulation year
    d13C_BGveg = "numeric", # The d13C value of roots vegetation in the last simulation year
    d13C_rhizodeposits = "numeric", # The d13C value of rhizodeposits in the last simulation year
    d13C_atmosphere = "numeric", # The d13C of atmospheric CO2 for every simulated year (calculated by the model based on data in the .csv file provided by the user)
    atm_ratio_13C_12C = "numeric", # The ratio of 13CO2 to 12CO2 in every simulated year (calculated by the model)
    CO2conc_atmosphere = "numeric", # The CO2 concentration in the atmosphere for every simulated year (calculated by the model based on data in the .csv file provided by the user)
    atm_ratio_14C_12C = "numeric", # The ratio of 14CO2 to 12CO2 in every simulated year (calculated by the model)
    d14C_atmosphere = "numeric", # The d14C of atmospheric CO2 for every simulated year (calculated by the model based on data in the .csv file provided by the user)
    d14C_atmosphere_years = "numeric", # The years for which d14CO2 data is provided (calculated by the model)
    AGveg_ratio_14C_12C = "numeric", # The ratio of 14C to 12C of aboveground vegetation for every simulated year (calculated by the model)
    ratioAtmosphic13CIsotopes = "numeric", # The ratio of 13CO2 to 12CO2 in every simulated year (calculated by the model)
    init_rhizo_13C_12C = "numeric", # Necessary for the calculation of the 13CO2 soil profile (the model choses this value based on other inputs)
    init_rhizo_14C_12C = "numeric", # Necessary for the calculation of the 14CO2 soil profile (the model choses this value based on other inputs)
    atmCO2_12C_MolPerCubicMeter = "numeric", # Calculated by the model
    CO2_12_init = "numeric", # Calculated by the model
    CO2_13_init = "numeric", # Calculated by the model
    CO2_14_init = "numeric" # Calculated by the model
    
  ))
}

# -----------------------------------------
# A function to create the vegetation class
# -----------------------------------------

createVegClass <- function(){
  
  setClass("vegObject", slots=list(
    
    i_agv = "numeric", # Annual C inputs from aboveground vegetation [kg C m-2 yr-1]
    i_bgveg_tot = "numeric", # Annual C inputs from belowground vegetation (roots and exudates combined) [kg C m-2 yr-1 down to 1m depth]
    
    calculateRootProfileExponentially = "numeric", # If 1, the belowground C inputs decrease exponentially with depth
    calculateRootProfileAssymptotically = "numeric", # If 1, the belowground C inputs decrease assymptotically with depth (the standard situation)
    
    e_depth_bgv = "numeric", # The e-folding depth for belowground C inputs (if calculateRootProfileExponentially == 1)
    betaRoots = "numeric", # The beta value for belowground C inputs (if calculateRootProfileAssymptotically == 1)
    
    AGveg_d13C = "numeric", # The d13C value of aboveground vegetation (= C inputs) for the last simulation year
    BGveg_d13C = "numeric", # The d13C value of roots (= C inputs) for the last simulation year
    d13C_rhizodeposits = "numeric",  # The d13C value of rhizodeposits (= C inputs) for the last simulation year
    
    fineRootBeta = "numeric", # The beta value for the depth profile of fine roots
    fineRootLength = "numeric", # Total length of fine roots down to 1 m depth [km m-3]
    exudationDist = "numeric", # Maximum distance of root exudation [mm]
    exudationK = "numeric", # Rate at which root exudate distance decreases with increasing root diameter [mm-1]
    
    i_bgveg_stage1 = "numeric", # In case 2 vegetation stages are run
    i_bgveg_stage2 = "numeric", # In case 2 vegetation stages are run
    
    annualRootCO2_stage1 = "numeric", # In case 2 vegetation stages are run
    annualRootCO2_stage2 = "numeric", # In case 2 vegetation stages are run
    
    rhizosphereVolumeFraction_vegStage1 = "numeric", # In case 2 vegetation stages are run
    rhizosphereVolumeFraction_vegStage2 = "numeric" # In case 2 vegetation stages are run
    
  ))
}