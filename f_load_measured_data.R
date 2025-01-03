# ------------------------
# Measured data are loaded
# ------------------------

loadMeasuredData <- function(site){
  
  if(site == "Hainich"){
  
    # The Excel file with the measured data
    # Every column contains some data, see data.frame below
    fileName <- "measuredData.xlsx"
    # The number of rows with soil data
    maxRow = 8
    
    # The data for the litter layer is loaded
    measuredData_litter <- data.frame(
      "Litter_Cstock" =  read_excel(fileName,"Hainich", "B2",col_names = "Litter_Cstock")[[1]],
      "Litter_d13C" =  read_excel(fileName,"Hainich", "D2",col_names = "Litter_Cstock")[[1]],
      "Litter_d14C" =  read_excel(fileName,"Hainich", "F2",col_names = "Litter_Cstock")[[1]]
    )
    
    # The data for the soil is loaded
    measuredData_soil <- data.frame(
      "Depth"         =  read_excel(fileName,"Hainich",paste("A3:A", as.character(maxRow), sep = ""),col_names = "Depth")[[1]],
      
      "fPOC_conc"     =  read_excel(fileName,"Hainich",paste("B3:B", as.character(maxRow), sep = ""),col_names = "Depth")[[1]],
      "MAOC_conc"     =  read_excel(fileName,"Hainich",paste("C3:C", as.character(maxRow), sep = ""),col_names = "Depth")[[1]],
      
      "fPOC_d13C"     =  read_excel(fileName,"Hainich",paste("D3:D", as.character(maxRow), sep = ""),col_names = "Depth")[[1]],
      "MAOC_d13C"     =  read_excel(fileName,"Hainich",paste("E3:E", as.character(maxRow), sep = ""),col_names = "Depth")[[1]],
      
      "fPOC_d14C"     =  read_excel(fileName,"Hainich",paste("F3:F", as.character(maxRow), sep = ""),col_names = "Depth")[[1]],
      "MAOC_d14C"     =  read_excel(fileName,"Hainich",paste("G3:G", as.character(maxRow), sep = ""),col_names = "Depth")[[1]],
      
      "lowerBound"     =  read_excel(fileName,"Hainich",paste("H3:H", as.character(maxRow), sep = ""),col_names = "Depth")[[1]],
      "upperBound"     =  read_excel(fileName,"Hainich",paste("I3:I", as.character(maxRow), sep = ""),col_names = "Depth")[[1]],
      
      "total_d13C"     =  read_excel(fileName,"Hainich",paste("M3:M", as.character(maxRow), sep = ""),col_names = "Depth")[[1]],
      "total_d14C"     =  read_excel(fileName,"Hainich",paste("N3:N", as.character(maxRow), sep = ""),col_names = "Depth")[[1]],
      
      "fPOC_conc_stDev"     =  read_excel(fileName,"Hainich",paste("T3:T", as.character(maxRow), sep = ""),col_names = "Depth")[[1]],
      "MAOC_conc_stDev"     =  read_excel(fileName,"Hainich",paste("U3:U", as.character(maxRow), sep = ""),col_names = "Depth")[[1]],
      
      "fPOC_d13C_stDev"     =  read_excel(fileName,"Hainich",paste("V3:V", as.character(maxRow), sep = ""),col_names = "Depth")[[1]],
      "MAOC_d13C_stDev"     =  read_excel(fileName,"Hainich",paste("W3:W", as.character(maxRow), sep = ""),col_names = "Depth")[[1]],
      
      "fPOC_d14C_stDev"     =  read_excel(fileName,"Hainich",paste("X3:X", as.character(maxRow), sep = ""),col_names = "Depth")[[1]],
      "MAOC_d14C_stDev"     =  read_excel(fileName,"Hainich",paste("Y3:Y", as.character(maxRow), sep = ""),col_names = "Depth")[[1]]
    )
    
    # -----------------------------------------
    # The "measurements" of Cmic are calculated
    # -----------------------------------------
    
    # The total OC%
    Ctot <- rowSums(cbind(measuredData_soil$fPOC_conc, measuredData_soil$MAOC_conc))
    
    # Rhizosphere: Cmic is 2 % of total OC
    Cmic_rhizo_conc <- Ctot*0.02
    
    # Bulk soil: Cmic is 1 % of total OC
    Cmic_bulk_conc <- Ctot*0.01
    
    # This is added to the data.frame
    measuredData_soil <- as.data.frame(c(measuredData_soil,
                           "Cmic_rhizo_conc" = list(Cmic_rhizo_conc),
                           "Cmic_bulk_conc" = list(Cmic_bulk_conc)))
    
    
    # ---------------------------------------------------------
    # Only measurements above the evaluation depth are retained
    # ---------------------------------------------------------
    
    # The maximum simulated depth
    maxDepth <- par@evalDepth
    
    # The row numbers with a depth above the maximum simulated depth
    rowNum <- which(measuredData_soil$Depth <= maxDepth)
    
    # Only these rows are retained
    measuredData_soil <- measuredData_soil[rowNum,]
    
  
  return(list("measuredData_litter" = measuredData_litter, "measuredData_soil" = measuredData_soil))

}
}