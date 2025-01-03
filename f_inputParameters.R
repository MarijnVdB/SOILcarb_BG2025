# In this script, the input parameters for the SOILcarb model are defined
# For info on parameters, see the file 'f_createClasses.R'

inputParameters <- function(par){
  
  if(par@site == "Hainich"){
    
      # -------------------------------
      # Parameters for all compartments
      # -------------------------------

      par@CUE_L = 0.3
      par@CUE_R = 0.3
      par@CUE_B = 0.3
      par@Db0_stage1 = 5.682698e-06
      par@Db_eFold_depth_stage1 = 0.15
      par@denityMin = data.frame(
                      "upperDepth" = c(0,5,10,20,30,40,50,70,90),
                      "lowerDepth" = c(5,10,20,30,40,50,60,80,100),
                      "bulkDensity" = c(0.79,0.99,1.20,1.31,1.41,1.48,1.56, 1.65, 1.7))
      par@advectionRate_polyC_stage1 = 0.80 # Constant, as no effect in sensitivity analysis
      par@fractIncorpCO2_litter = 0.011 # Akinyede et al. 2020 SBB
      par@fractIncorpCO2_soil = 0.011
      par@solubleMicFraction = 0.5
      par@beta = 2

      par@Db0_stage2 = 2.6e-5
      par@Db_eFold_depth_stage2 = 0.285
      
      par@CO2conc_fractionationFactor = 0.014

      # -------------------------------
      # Parameters for the litter layer
      # -------------------------------

      par@Vmax_L = 95
      par@KmB_POC_L_stage1 = 7.15
      par@KmB_DOC_L_stage1 = 3.53
      par@L_fractionLeachable = 0.4
      par@leachFract_stage1 = 0.20
      par@bioturbFract_stage1 = 0.10
      par@K_mic_L = 0.05  # Carrying capacity %

      # ------------------------------
      # Parameters for the rhizosphere
      # ------------------------------

      par@Fdoc_rhizo = 0.2857# Pausch and Kuzyakov (2010)
      
      par@VmaxD_root_R_stage1 = 9.980782e-01
      par@KmB_root_R_stage1 = 0.02 # MicToPoc_R
      
      par@VmaxU_BioAv_stage1 = 0.8
      par@KmU_Bioav_stage1 = 0.02 # MicToBioAv_R
      
      par@K_mic_R = .1 # Carrying capacity %
      
      par@bioToPoly <- 0.175

      # -------------------------------
      # Parameters for the mineral soil
      # -------------------------------

      par@VmaxD_M_stage1 = 1.124894e+02
      par@Km_depol_M_stage1 = 4.663196e-01 # !!! here defined as a % of total soil mass, equivalent to %SOC
      par@Vmax_ads = 9.891532e+02
      par@Km_ads = 1.550768e-01 # !!! here defined as a % of total soil mass, equivalent to %SOC
      par@kDes_init = 1.160612e-01
      par@CminMax = 0.083 # Relation from Georgiou 2022 (s+c% of 97 %)
      par@K_mic_B = 0.03 # Carrying capacity %

      # -------------------------------------------------------
      # Parameters for the calculation of the CO2 depth profile
      # -------------------------------------------------------

      par@Tair = 8
      par@airFilledPorosity_surface = 0.2
      par@atmPressure = 101.325
      par@dailyRootCO2Production = 0.00120

      # ------------------
      # Vegetation stage 1
      # ------------------
      
      vegData_stage1 <- new("vegObject")
      
      vegData_stage1@i_agv = 0.209 # Kutsch (2018)
      vegData_stage1@i_bgveg_tot = 0.324 # !!! Inputs down to 1 m depth !!!

      vegData_stage1@calculateRootProfileExponentially = 0
      vegData_stage1@calculateRootProfileAssymptotically = 1

      vegData_stage1@e_depth_bgv = 0
      vegData_stage1@betaRoots = 8.729074e-01

      vegData_stage1@AGveg_d13C = -29.3
      vegData_stage1@BGveg_d13C = -27.8
      vegData_stage1@d13C_rhizodeposits = -28.9

      vegData_stage1@fineRootBeta = vegData_stage1@betaRoots
      vegData_stage1@fineRootLength = 5.4
      vegData_stage1@exudationDist = 2
      vegData_stage1@exudationK = 1.5
      
      # ------------------
      # Vegetation stage 2
      # ------------------
      
      vegData_stage2 <- new("vegObject")
      
      vegData_stage2@i_agv = 0.232
      vegData_stage2@i_bgveg_tot = 0.390
      
      vegData_stage2@calculateRootProfileExponentially = 0
      vegData_stage2@calculateRootProfileAssymptotically = 1
      
      vegData_stage2@e_depth_bgv = 0
      vegData_stage2@betaRoots = 0.94
      
      vegData_stage2@AGveg_d13C = -30.2
      vegData_stage2@BGveg_d13C = -28.0
      
      vegData_stage2@fineRootBeta = 0.94
      vegData_stage2@fineRootLength = 5.4
      vegData_stage2@exudationDist = 2
      vegData_stage2@exudationK = 1.5

  return(list("par" = par,
              "vegData_stage1" = vegData_stage1,
              "vegData_stage2" = vegData_stage2))
  
  } # End if
} # End function
