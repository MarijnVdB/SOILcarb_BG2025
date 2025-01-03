# This script is used to extract the unique parameter combinations from the DE output

# Remove all variables
rm(list=ls(all = TRUE))
# Clear the console window
cat("\014")

library("mgcv")
library("ggplot2")

# The output file is loaded: chose the correct one
folder <- "DE output/C/"
# folder <- "DE output/C_d13C/"
# folder <- "DE output/C_d14C/"
# folder <- "DE output/C_d13C_d14C/"

filename = paste(folder, "DE_out.RData", sep = "")
load(file = filename)
rm(filename)

# The number of iterations is obtained
n <- length(outDE$member$storepop)

# The error is plotted in function of iteration number
df_error <- data.frame("iteration" = 1:(n+1),
                    "error" = outDE$member$bestvalit)

ggplot(df_error, aes(x = iteration, y = error)) +
  geom_line()

# All parameter combinations are retrieved, after the lowest error has been obtained
# paramCombi <- NA

nStart <- 1
for(ii in nStart:n){

  if(ii == nStart){
    paramCombi <- outDE$member$storepop[[ii]]
  } else{
    paramCombi <- rbind(paramCombi,outDE$member$storepop[[ii]])
  }
}

# Only unique rows with parameter values are retained
uniqueParam <- uniquecombs(paramCombi)

# These parameter values are saved: chose the correct folder
save(uniqueParam, file = "DE output/C/Unique_parameter_combinations/uniqueParam.RData")
# save(uniqueParam, file = "DE output/C_d13C/Unique_parameter_combinations/uniqueParam.RData")
# save(uniqueParam, file = "DE output/C_d14C/Unique_parameter_combinations/uniqueParam.RData")
# save(uniqueParam, file = "DE output/C_d13C_d14C/Unique_parameter_combinations/uniqueParam.RData")
