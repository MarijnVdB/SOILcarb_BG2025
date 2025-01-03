plotting <- function(res, d13C, d14C, par){
  
  # This function generates the plot objects, which are exported and plotted
  # in the main script
  
  # ----------------------------------------------------------------------------
  # The litter results are plotted
  # ----------------------------------------------------------------------------
  
  # A data frame with the data to plot
  df_litter <- data.frame("time" = seq(par@lastRealSimulationYear - par@numberOfSimulationYears+1, par@lastRealSimulationYear),
                          "POC_L" = res$POC_L,
                          "DOC_L" = res$DOC_L,
                          "MIC_L" = res$MIC_L,
                          "TOT_L" = res$C_litter,
                          "d13C_L" = d13C$C_L,
                          "d14C_L" = d14C$C_L)
  
  # -------------
  # Litter carbon
  # -------------
  
  pl1 <- ggplot() +
    geom_line(data = df_litter, aes(x = time, y = POC_L, color = "POC")) +
    geom_line(data = df_litter, aes(x = time, y = DOC_L, color = "DOC")) +
    geom_line(data = df_litter, aes(x = time, y = MIC_L, color = "Microbes")) +
    geom_line(data = df_litter, aes(x = time, y = TOT_L, color = "Total C")) +
    geom_point(data = measuredData_litter, aes(x = par@lastRealSimulationYear, y = Litter_Cstock, color = "Total C"), size = 2) +
    geom_point(data = measuredData_litter, aes(x = par@lastRealSimulationYear, y = Litter_Cstock * 0.33, color = "DOC"), size = 2) +
    geom_point(data = measuredData_litter, aes(x = par@lastRealSimulationYear, y = Litter_Cstock * 0.66, color = "POC"), size = 2) +
    geom_point(data = measuredData_litter, aes(x = par@lastRealSimulationYear, y = Litter_Cstock * 0.01, color = "Microbes"), size = 2) +
    labs(title="Litter organic carbon",
         y=expression(paste("OC (kg C ", m^-2, ")")),
         x="Year") +
    theme_classic() +
    scale_color_manual(name = NULL, values = c("POC" = "#cc4c02" ,
                                               "DOC" = "#0570b0",
                                               "Microbes" = "#88419d",
                                               "Total C" = "Black")) +
    scale_x_continuous(expand = c(0,0), limits = c(1850,par@lastRealSimulationYear+2)) +
    ylim(0,0.6)#0.35)
  
  # -----------
  # Litter d13C
  # -------------
  
  pl2 <- ggplot() +
    geom_line(data = df_litter, aes(x = time, y = d13C_L)) +
    geom_point(data = measuredData_litter, aes(x = par@lastRealSimulationYear, y = Litter_d13C), size = 2) +
    scale_x_continuous(expand = c(0,0), limits = c(1850,par@lastRealSimulationYear+2)) +
    labs(title=expression(paste("Litter ", delta^{13} , "C")),
         y=expression(paste(delta^{13}, "C (\u2030)")),
         x="Year") +
    theme_classic()
  
  # -------------
  # Litter d14C
  # -------------
  
  pl3 <- ggplot() +
    geom_line(data = df_litter, aes(x = time, y = d14C_L)) +
    geom_point(data = measuredData_litter, aes(x = par@lastRealSimulationYear, y = Litter_d14C), size = 2) +
    scale_x_continuous(expand = c(0,0), limits = c(1850,par@lastRealSimulationYear+2)) +
    labs(title=expression(paste("Litter ", Delta^{14}, "C")),
         y=expression(paste(Delta^{14}, "C (\u2030)")),
         x="Year") +
    theme_classic()
  
  # g <- ggarrange(pl1, pl2, pl3, nrow = 3, labels = c("(A)", "(B)", "(C)"), font.label=list(color="black",size=12))
  # ggsave("Figures for manuscript/Supplement_Figure_litter.png", g, device = "png", dpi = 300, units = "cm", width = 20, height = 12)
  
  # ----------------------------------------------------------------------------
  # The soil carbon results are plotted - all pools
  # ----------------------------------------------------------------------------
  
  # Colors are defined
  c_Cmin <- "#8c6bb1"
  c_POC <- "#007f73"
  c_BioAv <- "#99d8c9"
  c_DOC <- "#74a9cf"
  c_mic_rhizo <- "#fdbb84"
  c_mic_bulk <- "#fe9929"
  c_CO2 <- "#dd1c77"
  c_bulk <- "#a84f02"
  c_Ctot <- "black"
  
  # Constructing the data frame
  df_soil <- data.frame("depth" = par@midDepth_layer,
                        
                        "POC" = (res$POC_R[par@numberOfSimulationYears,] / par@soilMassArray) * 100,
                        "BIOAV" = (res$BIOAV_R[par@numberOfSimulationYears,] / par@soilMassArray) * 100,
                        "MIC_R" = (res$MIC_R[par@numberOfSimulationYears,] / par@soilMassArray) * 100,
                        "DOC" = (res$POLY_B[par@numberOfSimulationYears,] / par@soilMassArray) * 100,
                        "MIN" = (res$MIN_B[par@numberOfSimulationYears,] / par@soilMassArray) * 100,
                        "MIC_B" = (res$MIC_B[par@numberOfSimulationYears,] / par@soilMassArray) * 100,
                        "Cbulk" = (res$C_bulk[par@numberOfSimulationYears,] / par@soilMassArray) * 100,
                        
                        "d13C_POC" = d13C$POC_R[par@numberOfSimulationYears,],
                        "d13C_BIOAV" = d13C$BIOAV_R[par@numberOfSimulationYears,],
                        "d13C_MIC_R" = d13C$MIC_R[par@numberOfSimulationYears,],
                        "d13C_DOC" = d13C$POLY_B[par@numberOfSimulationYears,],
                        "d13C_MIN" = d13C$MIN_B[par@numberOfSimulationYears,],
                        "d13C_MIC_B" = d13C$MIC_B[par@numberOfSimulationYears,],
                        "d13C_Ctot" = d13C$soilC[par@numberOfSimulationYears,],
                        "d13C_Cbulk" = d13C$C_bulk[par@numberOfSimulationYears,],
                        "d13C_C02" = d13C$CO2_tot[par@numberOfSimulationYears,],
                        "d13C_R" = d13C$CO2_R[par@numberOfSimulationYears,],
                        "d13C_B" = d13C$CO2_B[par@numberOfSimulationYears,],
                        
                        "d14C_POC" = d14C$POC_R[par@numberOfSimulationYears,],
                        "d14C_BIOAV" = d14C$BIOAV_R[par@numberOfSimulationYears,],
                        "d14C_MIC_R" = d14C$MIC_R[par@numberOfSimulationYears,],
                        "d14C_DOC" = d14C$POLY_B[par@numberOfSimulationYears,],
                        "d14C_MIN" = d14C$MIN_B[par@numberOfSimulationYears,],
                        "d14C_MIC_B" = d14C$MIC_B[par@numberOfSimulationYears,],
                        "d14C_Ctot" = d14C$soilC[par@numberOfSimulationYears,],
                        "d14C_Cbulk" = d14C$C_bulk[par@numberOfSimulationYears,],
                        "d14C_C02" = d14C$CO2_tot[par@numberOfSimulationYears,],
                        
                        "TT_rhizo_day" = res$TT_rhizo_day,
                        "TT_bulk_day" = res$TT_bulk_day
                        )
  
  df_soil <- cbind(df_soil,
               "Ctot" = df_soil$POC + df_soil$BIOAV + df_soil$MIC_R + df_soil$DOC + 
                 df_soil$MIN + df_soil$MIC_B)
  
  lineWidth_meas <- 2
  lineWidth_other <- 1
  scatterSize<- 5
  fontSize<- 16

  # ----------------
  # OC concentration
  # ----------------
  
  pc1 <- ggplot() + 
    
    # Modeled data
    geom_path(data = df_soil, aes(x = BIOAV, y = depth, color = "BIOAV"), size = lineWidth_other) +
    geom_path(data = df_soil, aes(x = MIC_R, y = depth, color = "MIC_R"), size = lineWidth_other) +
    geom_path(data = df_soil, aes(x = DOC, y = depth, color = "DOC"), size = lineWidth_other) +
    geom_path(data = df_soil, aes(x = MIC_B, y = depth, color = "MIC_B"), size = lineWidth_other) +
    geom_path(data = df_soil, aes(x = MIN, y = depth, color = "MIN"), size = lineWidth_other) +
    geom_path(data = df_soil, aes(x = POC, y = depth, color = "POC"), size = lineWidth_meas) +
    geom_path(data = df_soil, aes(x = Cbulk, y = depth, color = "C_bulk"), size = lineWidth_meas) +
    geom_path(data = df_soil, aes(x = Ctot, y = depth, color = "Ctot"), size = lineWidth_meas-1) +
    
    # Measured data
    geom_point(data = measuredData_soil, aes(x = fPOC_conc, y = Depth), fill = c_POC, color = "black", size = scatterSize, pch = 21) +
    geom_point(data = measuredData_soil, aes(x = MAOC_conc, y = Depth), fill = c_bulk, color = "black", size = scatterSize, pch = 21) +
    geom_point(data = measuredData_soil, aes(x = fPOC_conc + MAOC_conc, y = Depth), size = scatterSize-2, color = c_Ctot) +
    geom_errorbarh(data = measuredData_soil, aes(y = Depth, xmin = fPOC_conc - fPOC_conc_stDev, xmax = fPOC_conc + fPOC_conc_stDev, height = .01), color = "black") +
    geom_errorbarh(data = measuredData_soil, aes(y = Depth, xmin = MAOC_conc - MAOC_conc_stDev, xmax = MAOC_conc + MAOC_conc_stDev, height = .01), color = "black") +
    
    scale_y_continuous(trans = "reverse", limits = c(.65, 0), expand = c(0,0), sec.axis = sec_axis(~ ., labels = NULL)) + 
    scale_x_continuous(position = "top", limits = c(0, 9), expand = c(0,0), sec.axis = sec_axis(~ ., labels = NULL)) + 
    labs(y="Depth (m)",
        x="Soil organic carbon (%)") +
    theme_classic() +
    scale_color_manual(name = NULL, 
                       values = c("POC" = c_POC,
                                 "BIOAV" = c_BioAv,
                                 "MIC_R" = c_mic_rhizo,
                                 "DOC" = c_DOC,
                                 "MIN" = c_Cmin,
                                 "MIC_B" = c_mic_bulk,
                                 "C_bulk" = c_bulk,
                                 "Ctot" = c_Ctot),
                       labels = c("POC" = expression("C"["POC-r"]),
                                  "BIOAV" = expression("C"["bioav-r"]),
                                  "MIC_R" = expression("C"["mic-r"]),
                                  "DOC" = expression("C"["DOC-b"]),
                                  "MIN" = expression("C"["min-b"]),
                                  "MIC_B" = expression("C"["MIC-b"]),
                                  "C_bulk" = "Bulk soil C",
                                  "Ctot" = "Total carbon")) +
    scale_fill_manual(name = NULL, values = c("POC" = c_POC,
                                               "BIOAV" = c_BioAv,
                                               "MIC_R" = c_mic_rhizo,
                                               "DOC" = c_DOC,
                                               "MIN" = c_Cmin,
                                               "MIC_B" = c_mic_bulk,
                                               "C_bulk" = c_bulk,
                                               "Ctot" = "white")) +
    theme(legend.text.align = 0,
          legend.text=element_text(size = fontSize),
          legend.box.background = element_rect(colour = "black", fill = NULL, linetype = "solid"),
          # legend.box.margin = margin(1,2,1,1,"mm"),
          legend.position = c(0.68,0.30),
          axis.text=element_text(size = fontSize),
          axis.title=element_text(size = fontSize),
          plot.margin = margin(0,10,0,10, "pt")) #top, right, bottom, left
  
  # ----
  # d13C
  # ----
  
  pc2 <- ggplot() + 
    
    # Modeled data
    geom_path(data = df_soil, aes(x = d13C_BIOAV, y = depth, color = "BIOAV"), size = lineWidth_other) +
    geom_path(data = df_soil, aes(x = d13C_MIC_R, y = depth, color = "MIC_R"), size = lineWidth_other) +
    geom_path(data = df_soil, aes(x = d13C_DOC, y = depth, color = "DOC"), size = lineWidth_other) +
    geom_path(data = df_soil, aes(x = d13C_MIN, y = depth, color = "MIN"), size = lineWidth_other) +
    geom_path(data = df_soil, aes(x = d13C_MIC_B, y = depth, color = "MIC_B"), size = lineWidth_other) +
    geom_path(data = df_soil, aes(x = d13C_POC, y = depth, color = "POC"), size = lineWidth_meas) +
    geom_path(data = df_soil, aes(x = d13C_Cbulk, y = depth, color = "C_bulk"), size = lineWidth_meas) +
    geom_path(data = df_soil, aes(x = d13C_Ctot, y = depth, color = "Ctot"), size = lineWidth_meas-1) +
    geom_path(data = df_soil, aes(x = d13C_C02, y = depth, color = "C_CO2"), size = lineWidth_other, linetype = "dashed") +
    
    # Measured data
    geom_errorbarh(data = measuredData_soil, aes(y = Depth, xmin = fPOC_d13C - fPOC_d13C_stDev, xmax = fPOC_d13C + fPOC_d13C_stDev, height = .01), color = "black") +
    geom_errorbarh(data = measuredData_soil, aes(y = Depth, xmin = MAOC_d13C - MAOC_d13C_stDev, xmax = MAOC_d13C + MAOC_d13C_stDev, height = .01), color = "black") +
    geom_point(data = measuredData_soil, aes(x = fPOC_d13C, y = Depth), fill = c_POC, color = "black", size = scatterSize, pch = 21) +
    geom_point(data = measuredData_soil, aes(x = MAOC_d13C, y = Depth), fill = c_bulk, color = "black", size = scatterSize, pch = 21) +
    geom_point(data = measuredData_soil, aes(x = total_d13C, y = Depth), color = c_Ctot, size = scatterSize-2) +
    
    scale_y_continuous(trans = "reverse", limits = c(.65, 0), expand = c(0,0), sec.axis = sec_axis(~ ., labels = NULL)) + 
    scale_x_continuous(position = "top", limits = c(-29, -24), expand = c(0,0), sec.axis = sec_axis(~ ., labels = NULL)) + 
    labs(y="Depth (m)",
         x=expression(paste(delta^{13}, "C (\u2030)"))) +
    theme_classic() +
    scale_color_manual(name = NULL, 
                       values = c("POC" = c_POC,
                                               "BIOAV" = c_BioAv,
                                               "MIC_R" = c_mic_rhizo,
                                               "DOC" = c_DOC,
                                               "MIN" = c_Cmin,
                                               "MIC_B" = c_mic_bulk,
                                               "C_bulk" = c_bulk,
                                               "Ctot" = c_Ctot,
                                               "C_CO2" = c_CO2),
                       labels = c("POC" = expression("C"["POC-r"]),
                                  "BIOAV" = expression("C"["bioav-r"]),
                                  "MIC_R" = expression("C"["mic-r"]),
                                  "DOC" = expression("C"["DOC-b"]),
                                  "MIN" = expression("C"["min-b"]),
                                  "MIC_B" = expression("C"["mic-b"]),
                                  "C_bulk" = "Bulk soil C",
                                  "Ctot" = "Total carbon",
                                  "C_CO2" = "CO2")) +
    scale_fill_manual(name = NULL, values = c("POC" = c_POC,
                                              "BIOAV" = c_BioAv,
                                              "MIC_R" = c_mic_rhizo,
                                              "DOC" = c_DOC,
                                              "MIN" = c_Cmin,
                                              "MIC_B" = c_mic_bulk,
                                              "C_bulk" = c_bulk,
                                              "Ctot" = "white")) +
    theme(legend.position = "none",
          axis.text=element_text(size = fontSize),
          axis.title=element_text(size = fontSize),
          plot.margin = margin(0,10,0,10, "pt")) #top, right, bottom, left
  
  pc3 <- ggplot() + 
    
    # Modeled data
    geom_path(data = df_soil, aes(x = d14C_BIOAV, y = depth, color = "BIOAV"), size = lineWidth_other) +
    geom_path(data = df_soil, aes(x = d14C_MIC_R, y = depth, color = "MIC_R"), size = lineWidth_other) +
    geom_path(data = df_soil, aes(x = d14C_DOC, y = depth, color = "DOC"), size = lineWidth_other) +
    geom_path(data = df_soil, aes(x = d14C_MIN, y = depth, color = "MIN"), size = lineWidth_other) +
    geom_path(data = df_soil, aes(x = d14C_MIC_B, y = depth, color = "MIC_B"), size = lineWidth_other) +
    geom_path(data = df_soil, aes(x = d14C_POC, y = depth, color = "POC"), size = lineWidth_meas) +
    geom_path(data = df_soil, aes(x = d14C_Cbulk, y = depth, color = "C_bulk"), size = lineWidth_meas) +
    # geom_path(data = df_soil, aes(x = d14C_C02, y = depth, color = "C_CO2"), size = lineWidth_other, linetype = "dashed") +
    geom_path(data = df_soil, aes(x = d14C_Ctot, y = depth, color = "Ctot"), size = lineWidth_meas-1) +
    
    # Measured data
    geom_errorbarh(data = measuredData_soil, aes(y = Depth, xmin = fPOC_d14C - fPOC_d14C_stDev, xmax = fPOC_d14C + fPOC_d14C_stDev, height = .01), color = "black") +
    geom_errorbarh(data = measuredData_soil, aes(y = Depth, xmin = MAOC_d14C - MAOC_d14C_stDev, xmax = MAOC_d14C + MAOC_d14C_stDev, height = .01), color = "black") +
    geom_point(data = measuredData_soil, aes(x = fPOC_d14C, y = Depth), fill = c_POC, color = "black", size = scatterSize, pch = 21) +
    geom_point(data = measuredData_soil, aes(x = MAOC_d14C, y = Depth), fill = c_bulk, color = "black", size = scatterSize, pch = 21) +
    geom_point(data = measuredData_soil, aes(x = total_d14C, y = Depth), color = c_Ctot, size = scatterSize-2) +
    
    scale_y_continuous(trans = "reverse", limits = c(.65, 0), expand = c(0,0), sec.axis = sec_axis(~ ., labels = NULL)) + 
    scale_x_continuous(position = "top", limits = c(-610, 210), expand = c(0,0), sec.axis = sec_axis(~ ., labels = NULL)) + 
    labs(y="Depth (m)",
         x=expression(paste(Delta^{14}, "C (\u2030)"))) +
    theme_classic() +
    scale_color_manual(name = NULL, 
                       values = c("POC" = c_POC,
                                 "BIOAV" = c_BioAv,
                                 "MIC_R" = c_mic_rhizo,
                                 "DOC" = c_DOC,
                                 "MIN" = c_Cmin,
                                 "MIC_B" = c_mic_bulk,
                                 "C_bulk" = c_bulk,
                                 "Ctot" = c_Ctot,
                                 "C_CO2" = c_CO2),
                       labels = c("POC" = expression("C"["POC-r"]),
                                  "BIOAV" = expression("C"["bioav-r"]),
                                  "MIC_R" = expression("C"["mic-r"]),
                                  "DOC" = expression("C"["DOC-b"]),
                                  "MIN" = expression("C"["min-b"]),
                                  "MIC_B" = expression("C"["mic-b"]),
                                  "C_bulk" = "Bulk soil C",
                                  "Ctot" = "Total carbon",
                                  "C_CO2" = "CO2")) +
    scale_fill_manual(name = NULL, values = c("POC" = c_POC,
                                              "BIOAV" = c_BioAv,
                                              "MIC_R" = c_mic_rhizo,
                                              "DOC" = c_DOC,
                                              "MIN" = c_Cmin,
                                              "MIC_B" = c_mic_bulk,
                                              "C_bulk" = c_bulk,
                                              "Ctot" = "white")) +
    theme(legend.position = "none",
          axis.text=element_text(size = fontSize),
          axis.title=element_text(size = fontSize),
          plot.margin = margin(0,15,0,10, "pt")) #top, right, bottom, left

  # A new plotting window
  # quartz(title = "Depth profiles", width = 16*0.39, height = 7*0.39)

  # The plot is saved
  # g <- ggarrange(pc1,pc2,pc3, nrow = 1, labels = c("(A)", "(B)", "(C)"), font.label=list(color="black",size=12))
  # ggsave("Figures for manuscript/Supplement_allModelledPools.png", g, device = "png", dpi = 300, units = "cm", width = 32, height = 12)
  
  # ----------------------------------------------------------------------------
  # The soil carbon results are plotted - measured pools (for manuscript)
  # ----------------------------------------------------------------------------
  
  # ----------------
  # OC concentration
  # ----------------
  
  p1 <- ggplot() +

    # Modeled data
    geom_path(data = df_soil, aes(x = POC, y = depth, color = "POC"), size = lineWidth_meas) +
    geom_path(data = df_soil, aes(x = Cbulk, y = depth, color = "C_bulk"), size = lineWidth_meas) +
    geom_path(data = df_soil, aes(x = Ctot, y = depth, color = "Ctot"), size = lineWidth_meas-1) +

    # Measured data
    geom_point(data = measuredData_soil, aes(x = fPOC_conc, y = Depth), fill = c_POC, color = "black", size = scatterSize, pch = 21) +
    geom_point(data = measuredData_soil, aes(x = MAOC_conc, y = Depth), fill = c_bulk, color = "black", size = scatterSize, pch = 21) +
    geom_errorbarh(data = measuredData_soil, aes(y = Depth, xmin = fPOC_conc - fPOC_conc_stDev, xmax = fPOC_conc + fPOC_conc_stDev, height = .01), color = c_POC) +
    geom_errorbarh(data = measuredData_soil, aes(y = Depth, xmin = MAOC_conc - MAOC_conc_stDev, xmax = MAOC_conc + MAOC_conc_stDev, height = .01), color = c_bulk) +
    geom_point(data = measuredData_soil, aes(x = fPOC_conc + MAOC_conc, y = Depth), size = scatterSize-2, color = c_Ctot) +

    scale_y_continuous(trans = "reverse", limits = c(.65, 0), expand = c(0,0), sec.axis = sec_axis(~ ., labels = NULL)) +
    scale_x_continuous(position = "top", limits = c(0, 9), expand = c(0,0), sec.axis = sec_axis(~ ., labels = NULL)) +
    labs(y="Depth (m)",
         x="Soil organic carbon (%)") +
    theme_classic() +
    scale_color_manual(name = NULL,
                       values = c("POC" = c_POC,
                                  "C_bulk" = c_bulk,
                                  "Ctot" = c_Ctot),
                       labels = c("POC" = expression("C"["POC-r"]),
                                  "C_bulk" = "Bulk soil C",
                                  "Ctot" = "Total carbon")) +
    scale_fill_manual(name = NULL, values = c("POC" = c_POC,
                                              "C_bulk" = c_bulk,
                                              "Ctot" = "white")) +
    theme(legend.text.align = 0,
          legend.text=element_text(size = 16),
          legend.box.background = element_rect(colour = "black", fill = NULL, linetype = "solid"),
          legend.position = c(0.67,0.15),
          axis.text=element_text(size = fontSize),
          axis.title=element_text(size = fontSize),
          plot.margin = margin(0,10,0,10, "pt")) #top, right, bottom, left
  
  # ----
  # d13C
  # ----
  
  p2 <- ggplot() +

    # Modeled data
    geom_path(data = df_soil, aes(x = d13C_POC, y = depth, color = "POC"), size = lineWidth_meas) +
    geom_path(data = df_soil, aes(x = d13C_Cbulk, y = depth, color = "C_bulk"), size = lineWidth_meas) +
    geom_path(data = df_soil, aes(x = d13C_Ctot, y = depth, color = "Ctot"), size = lineWidth_meas-1) +

    # Measured data
    geom_errorbarh(data = measuredData_soil, aes(y = Depth, xmin = fPOC_d13C - fPOC_d13C_stDev, xmax = fPOC_d13C + fPOC_d13C_stDev, height = .01), color = c_POC) +
    geom_errorbarh(data = measuredData_soil, aes(y = Depth, xmin = MAOC_d13C - MAOC_d13C_stDev, xmax = MAOC_d13C + MAOC_d13C_stDev, height = .01), color = c_bulk) +
    geom_point(data = measuredData_soil, aes(x = fPOC_d13C, y = Depth), fill = c_POC, color = "black", size = scatterSize, pch = 21) +
    geom_point(data = measuredData_soil, aes(x = MAOC_d13C, y = Depth), fill = c_bulk, color = "black", size = scatterSize, pch = 21) +
    geom_point(data = measuredData_soil, aes(x = total_d13C, y = Depth), color = c_Ctot, size = scatterSize-2) +

    scale_y_continuous(trans = "reverse", limits = c(.65, 0), expand = c(0,0), sec.axis = sec_axis(~ ., labels = NULL)) +
    scale_x_continuous(position = "top", limits = c(-29, -24), expand = c(0,0), sec.axis = sec_axis(~ ., labels = NULL)) +
    scale_color_manual(name = NULL,
                       values = c("POC" = c_POC,
                                  "C_bulk" = c_bulk,
                                  "Ctot" = c_Ctot),
                       labels = c("POC" = expression("C"["POC-r"]),
                                  "C_bulk" = "Bulk soil C",
                                  "Ctot" = "Total carbon")) +
    scale_fill_manual(name = NULL, values = c("POC" = c_POC,
                                              "C_bulk" = c_bulk,
                                              "Ctot" = "white")) +
    labs(y="Depth (m)",
         x=expression(paste(delta^{13}, "C (\u2030)"))) +
    theme_classic() +
    theme(axis.text=element_text(size = fontSize),
          axis.title=element_text(size = fontSize),
          plot.margin = margin(0,10,0,10, "pt"),
          legend.position = "none") #top, right, bottom, left
  
  p3 <- ggplot() +

    # Modeled data
    geom_path(data = df_soil, aes(x = d14C_POC, y = depth, color = "POC"), size = lineWidth_meas) +
    geom_path(data = df_soil, aes(x = d14C_Cbulk, y = depth, color = "C_bulk"), size = lineWidth_meas) +
    geom_path(data = df_soil, aes(x = d14C_Ctot, y = depth, color = "Ctot"), size = lineWidth_meas-1) +

    # Measured data
    geom_errorbarh(data = measuredData_soil, aes(y = Depth, xmin = fPOC_d14C - fPOC_d14C_stDev, xmax = fPOC_d14C + fPOC_d14C_stDev, height = .01), color = c_POC) +
    geom_errorbarh(data = measuredData_soil, aes(y = Depth, xmin = MAOC_d14C - MAOC_d14C_stDev, xmax = MAOC_d14C + MAOC_d14C_stDev, height = .01), color = c_bulk) +
    geom_point(data = measuredData_soil, aes(x = fPOC_d14C, y = Depth), fill = c_POC, color = "black", size = scatterSize, pch = 21) +
    geom_point(data = measuredData_soil, aes(x = MAOC_d14C, y = Depth), fill = c_bulk, color = "black", size = scatterSize, pch = 21) +
    geom_point(data = measuredData_soil, aes(x = total_d14C, y = Depth), color = c_Ctot, size = scatterSize-2) +

    scale_y_continuous(trans = "reverse", limits = c(.65, 0), expand = c(0,0), sec.axis = sec_axis(~ ., labels = NULL)) +
    scale_x_continuous(position = "top", limits = c(-610, 210), expand = c(0,0), sec.axis = sec_axis(~ ., labels = NULL)) +
    scale_color_manual(name = NULL,
                       values = c("POC" = c_POC,
                                  "C_bulk" = c_bulk,
                                  "Ctot" = c_Ctot),
                       labels = c("POC" = expression("C"["POC-r"]),
                                  "C_bulk" = "Bulk soil C",
                                  "Ctot" = "Total carbon")) +
    scale_fill_manual(name = NULL, values = c("POC" = c_POC,
                                              "C_bulk" = c_bulk,
                                              "Ctot" = "white")) +
    labs(y="Depth (m)",
         x=expression(paste(Delta^{14}, "C (\u2030)"))) +
    theme_classic() +
    theme(axis.text=element_text(size = fontSize),
          axis.title=element_text(size = fontSize),
          plot.margin = margin(0,15,0,10, "pt"),
          legend.position = "none") #top, right, bottom, left
  
  # A new plotting window
  # quartz(title = "Depth profiles", width = 16*0.39, height = 7*0.39)
  
  # The plot is saved
  g <- ggarrange(p1,p2,p3, nrow = 1, labels = c("(a)", "(b)", "(c)"), font.label=list(color="black",size=18, face = "plain"))
  ggsave("Figures for manuscript/C_d13C_d14C_depthProfiles_measuredPools.png", g, device = "png", dpi = 300, units = "cm", width = 32, height = 12)
  
  # g <- ggarrange(p1,p2,p3, nrow = 1, labels = c("(A)", "(B)", "(C)"), font.label=list(color="black",size=12))
  # ggsave("Figures for manuscript/Supplement_d13C_onlyAGvegMixing.png", g, device = "png", dpi = 300, units = "cm", width = 32, height = 12)
  
  # g <- ggarrange(p1,p2,p3, nrow = 1, labels = c("(A)", "(B)", "(C)"), font.label=list(color="black",size=12))
  # ggsave("Figures for manuscript/Supplement_d13C_AGvegMixing_and_d13CO2noS.png", g, device = "png", dpi = 300, units = "cm", width = 32, height = 12)
  
  # ----------------------------------------------------------------------------
  # The soil carbon results are plotted - CO2 isotopes - for manuscript
  # ----------------------------------------------------------------------------
  
  # A new plotting window
  # quartz(title = "Depth profiles", width = 16*0.39, height = 7*0.39)
  
  lineWidth_meas <- 2
  lineWidth_other <- 1
  scatterSize<- 5
  fontSize<- 16
  
  # ----
  # d13C
  # ----
  
  p1 <- ggplot() + 
    
    # Modeled data
    geom_path(data = df_soil, aes(x = d13C_POC, y = depth, color = "POC"), size = lineWidth_meas) +
    geom_path(data = df_soil, aes(x = d13C_Cbulk, y = depth, color = "C_bulk"), size = lineWidth_meas) +
    geom_path(data = df_soil, aes(x = d13C_Ctot, y = depth, color = "Ctot"), size = lineWidth_meas-1) +
    geom_path(data = df_soil, aes(x = d13C_C02, y = depth, color = "C_CO2"), size = lineWidth_other, linetype = "dashed") +
    
    # Measured data
    geom_errorbarh(data = measuredData_soil, aes(y = Depth, xmin = fPOC_d13C - fPOC_d13C_stDev, xmax = fPOC_d13C + fPOC_d13C_stDev, height = .01), color = "black") +
    geom_errorbarh(data = measuredData_soil, aes(y = Depth, xmin = MAOC_d13C - MAOC_d13C_stDev, xmax = MAOC_d13C + MAOC_d13C_stDev, height = .01), color = "black") +
    geom_point(data = measuredData_soil, aes(x = fPOC_d13C, y = Depth), color = c_POC, size = scatterSize) +
    geom_point(data = measuredData_soil, aes(x = MAOC_d13C, y = Depth), color = c_bulk, size = scatterSize) +
    geom_point(data = measuredData_soil, aes(x = total_d13C, y = Depth), color = c_Ctot, size = scatterSize, shape = 21) +
    
    scale_y_continuous(trans = "reverse", limits = c(.65, 0), expand = c(0,0), sec.axis = sec_axis(~ ., labels = NULL)) + 
    scale_x_continuous(position = "top", limits = c(-29, -22), expand = c(0,0), sec.axis = sec_axis(~ ., labels = NULL)) + 
    labs(y="Depth (m)",
         x=expression(paste(delta^{13}, "C (\u2030)"))) +
    theme_classic() +
    scale_color_manual(name = NULL, 
                       values = c("POC" = c_POC,
                                  "C_bulk" = c_bulk,
                                  "Ctot" = c_Ctot,
                                  "C_CO2" = c_CO2),
                       labels = c("POC" = expression("C"["POC-r"]),
                                  "C_bulk" = "Bulk soil C",
                                  "Ctot" = "Total carbon",
                                  "C_CO2" = "CO2")) +
    scale_fill_manual(name = NULL, values = c("POC" = c_POC,
                                              "C_bulk" = c_bulk,
                                              "Ctot" = c_Ctot,
                                              "C_CO2" = c_CO2)) +
    theme(legend.text.align = 0,
          legend.text=element_text(size = fontSize),
          legend.box.background = element_rect(colour = "black", fill = NULL, linetype = "solid"),
          legend.box.margin = margin(1,2,1,1,"mm"),
          axis.text=element_text(size = fontSize),
          axis.title=element_text(size = fontSize),
          legend.position ="none",
          plot.margin = margin(0,10,0,10, "pt")) #top, right, bottom, left
  
  p2 <- ggplot() + 
    
    # Modeled data
    geom_path(data = df_soil, aes(x = d14C_POC, y = depth, color = "POC"), size = lineWidth_meas) +
    geom_path(data = df_soil, aes(x = d14C_Cbulk, y = depth, color = "C_bulk"), size = lineWidth_meas) +
    geom_path(data = df_soil, aes(x = d14C_Ctot, y = depth, color = "Ctot"), size = lineWidth_meas-1) +
    geom_path(data = df_soil, aes(x = d14C_C02, y = depth, color = "C_CO2"), size = lineWidth_other, linetype = "dashed") +
    
    # Measured data
    geom_errorbarh(data = measuredData_soil, aes(y = Depth, xmin = fPOC_d14C - fPOC_d14C_stDev, xmax = fPOC_d14C + fPOC_d14C_stDev, height = .01), color = "black") +
    geom_errorbarh(data = measuredData_soil, aes(y = Depth, xmin = MAOC_d14C - MAOC_d14C_stDev, xmax = MAOC_d14C + MAOC_d14C_stDev, height = .01), color = "black") +
    geom_point(data = measuredData_soil, aes(x = fPOC_d14C, y = Depth), color = c_POC, size = scatterSize) +
    geom_point(data = measuredData_soil, aes(x = MAOC_d14C, y = Depth), color = c_bulk, size = scatterSize) +
    
    scale_y_continuous(trans = "reverse", limits = c(.65, 0), expand = c(0,0), sec.axis = sec_axis(~ ., labels = NULL)) + 
    scale_x_continuous(position = "top", limits = c(-610, 210), expand = c(0,0), sec.axis = sec_axis(~ ., labels = NULL)) + 
    labs(y="Depth (m)",
         x=expression(paste(Delta^{14}, "C (\u2030)"))) +
    theme_classic() +
    scale_color_manual(name = NULL, 
                       values = c("POC" = c_POC,
                                  "C_bulk" = c_bulk,
                                  "Ctot" = c_Ctot,
                                  "C_CO2" = c_CO2),
                       labels = c("POC" = expression("C"["POC-r"]),
                                  "C_bulk" = "Bulk soil C",
                                  "Ctot" = "Total carbon",
                                  "C_CO2" = "CO2")) +
    scale_fill_manual(name = NULL, values = c("POC" = c_POC,
                                              "C_bulk" = c_bulk,
                                              "Ctot" = c_Ctot,
                                              "C_CO2" = c_CO2)) +
    theme(legend.text.align = 0,
          legend.text=element_text(size = 12),
          legend.box.background = element_rect(colour = "black", fill = NULL, linetype = "solid"),
          legend.box.margin = margin(.5,.5,.5,.5,"mm"),
          axis.text=element_text(size = fontSize),
          axis.title=element_text(size = fontSize),
          plot.margin = margin(0,15,0,10, "pt"), #top, right, bottom, left
          plot.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white"),
          legend.position = c(0.3,0.8)) 
  
  # A new plotting window
  # quartz(title = "Depth profiles", width = 12*0.39, height = 7*0.39)
  
  # ggarrange(p1, p2, nrow = 1, common.legend = TRUE, legend = "right", labels = c("(A)", "(B)"), font.label=list(color="black",size=12))
  
  # The plot is saved
  # g <- ggarrange(p1, p2, nrow = 1, labels = c("(A)", "(B)"), font.label=list(color="black",size=12))
  # ggsave("Figures for manuscript/Supplement_Figure_CO2.png", g, device = "png", dpi = 300, units = "cm", width = 20, height = 12)
  
  # The plot objects are returned from the function
  return(list(pl1,pl2,pl3,pc1,pc2,pc3))
  
}
