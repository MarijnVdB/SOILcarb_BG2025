# This file is used to export the .RDATA file from the SAFE sensitivity run
# as .csv so that it can be imported in Matlab for further analyses

# Remove all variables
rm(list=ls(all = TRUE))
# Clear the console window
cat("\014")

folder <- "SAFE output/PAWN output - all parameters/"

# The file is loaded
out <- load(file = paste(folder, "output_SAFE.RData", sep = ""))

# The files is saved as .csv in the same folder
write.csv(output_SAFE, paste(folder, "/output_SAFE.csv", sep = ""))
