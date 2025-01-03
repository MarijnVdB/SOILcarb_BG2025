# This script is used to save the model output for the simple sensitivity analysis

# The output data is stored in a data.frame
df_soil <- data.frame("depth" = par@midDepth_layer,
                      
                      "POC" = (res$POC_R[par@numberOfSimulationYears,] / par@soilMassArray) * 100,
                      "BIOAV" = (res$BIOAV_R[par@numberOfSimulationYears,] / par@soilMassArray) * 100,
                      "MIC_R" = (res$MIC_R[par@numberOfSimulationYears,] / par@soilMassArray) * 100,
                      "POLY" = (res$POLY_B[par@numberOfSimulationYears,] / par@soilMassArray) * 100,
                      "MIN" = (res$MIN_B[par@numberOfSimulationYears,] / par@soilMassArray) * 100,
                      "MIC_B" = (res$MIC_B[par@numberOfSimulationYears,] / par@soilMassArray) * 100,
                      "Cbulk" = (res$C_bulk[par@numberOfSimulationYears,] / par@soilMassArray) * 100,
                      
                      "d13C_POC" = d13C$POC_R[par@numberOfSimulationYears,],
                      "d13C_BIOAV" = d13C$BIOAV_R[par@numberOfSimulationYears,],
                      "d13C_MIC_R" = d13C$MIC_R[par@numberOfSimulationYears,],
                      "d13C_POLY" = d13C$POLY_B[par@numberOfSimulationYears,],
                      "d13C_MIN" = d13C$MIN_B[par@numberOfSimulationYears,],
                      "d13C_MIC_B" = d13C$MIC_B[par@numberOfSimulationYears,],
                      "d13C_Ctot" = d13C$soilC[par@numberOfSimulationYears,],
                      "d13C_Cbulk" = d13C$C_bulk[par@numberOfSimulationYears,],
                      "d13C_C02" = d13C$CO2_tot[par@numberOfSimulationYears,],
                      
                      "d14C_POC" = d14C$POC_R[par@numberOfSimulationYears,],
                      "d14C_BIOAV" = d14C$BIOAV_R[par@numberOfSimulationYears,],
                      "d14C_MIC_R" = d14C$MIC_R[par@numberOfSimulationYears,],
                      "d14C_POLY" = d14C$POLY_B[par@numberOfSimulationYears,],
                      "d14C_MIN" = d14C$MIN_B[par@numberOfSimulationYears,],
                      "d14C_MIC_B" = d14C$MIC_B[par@numberOfSimulationYears,],
                      "d14C_Ctot" = d14C$soilC[par@numberOfSimulationYears,],
                      "d14C_Cbulk" = d14C$C_bulk[par@numberOfSimulationYears,],
                      "d14C_C02" = d14C$CO2_tot[par@numberOfSimulationYears,]
)

df_soil <- cbind(df_soil,
                 "Ctot" = df_soil$POC + df_soil$BIOAV + df_soil$MIC_R + df_soil$POLY + 
                   df_soil$MIN + df_soil$MIC_B)

# This data.frame is saved
folder <- "Simple parameter sensitivity/Data/"

# fileName <- "averageResult"

# fileName <- "d13C_AGveg_29pt9"
# fileName <- "d13C_AGveg_28pt9"
# 
# fileName <- "d13C_BGveg_28pt3"
# fileName <- "d13C_BGveg_27pt3"
# 
# fileName <- "d13C_rhizodeposits_29pt4"
# fileName <- "d13C_rhizodeposits_28pt4"
# 
# fileName <- "f_CO2_uptake_0"
# fileName <- "f_CO2_uptake_0pt05"
# 
# fileName <- "pCO2_fact_0pt0108"
fileName <- "pCO2_fact_0pt0172"

save(df_soil, file = paste(folder, fileName, ".RData", sep = ""))
