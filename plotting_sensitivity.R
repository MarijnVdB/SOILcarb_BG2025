# Remove all variables
rm(list=ls(all = TRUE))
# Clear the console window
cat("\014")
# Close all plot windows
graphics.off()

library("ggplot2")
library("gridExtra")
library("ggpubr")
library("dplyr")

# ------------------------------------------------------------------------------
# The data is loaded
# ------------------------------------------------------------------------------

# --------------------------
# Local sensitivity analysis
# --------------------------

folder <- "Simple parameter sensitivity/Data/"

# A function to load the data
loadData <- function(folder, fileName){
  load(paste(folder, fileName, ".Rdata", sep = ""))
  return(df_soil)
}

fileName <- "averageResult"
assign(fileName, loadData(folder, fileName))

fileName <- "d13C_AGveg_29pt9"
assign(fileName, loadData(folder, fileName))

fileName <- "d13C_AGveg_28pt9"
assign(fileName, loadData(folder, fileName))

fileName <- "d13C_BGveg_28pt3"
assign(fileName, loadData(folder, fileName))

fileName <- "d13C_BGveg_27pt3"
assign(fileName, loadData(folder, fileName))

fileName <- "d13C_rhizodeposits_29pt4"
assign(fileName, loadData(folder, fileName))

fileName <- "d13C_rhizodeposits_28pt4"
assign(fileName, loadData(folder, fileName))

fileName <- "f_CO2_uptake_0"
assign(fileName, loadData(folder, fileName))

fileName <- "f_CO2_uptake_0pt05"
assign(fileName, loadData(folder, fileName))

fileName <- "pCO2_fact_0pt0108"
assign(fileName, loadData(folder, fileName))

fileName <- "pCO2_fact_0pt0172"
assign(fileName, loadData(folder, fileName))

# ------------------------------------------------------------------------------
# The depth profiles are interpolated to a resolution of 1cm to make plotting
# easier
# ------------------------------------------------------------------------------

maxDepth <- floor(max(d13C_AGveg_28pt9$depth)*100)/100
depth_interp <- seq(0.005,maxDepth-0.005,0.01)

d13C_AGveg_28pt9_interp <- data.frame("depth" = depth_interp,
                                      "d13C_Ctot" = approx(d13C_AGveg_28pt9$depth,d13C_AGveg_28pt9$d13C_Ctot,depth_interp)$y)

# ggplot() +
#   geom_point(data = d13C_AGveg_28pt9_interp, aes(x = d13C_Ctot, y = depth)) +
#   geom_point(data = d13C_AGveg_28pt9, aes(x = d13C_Ctot, y = depth), color = "red")


d13C_AGveg_29pt9_interp <- data.frame("depth" = depth_interp,
                                      "d13C_Ctot" = approx(d13C_AGveg_29pt9$depth,d13C_AGveg_29pt9$d13C_Ctot,depth_interp)$y)

d13C_BGveg_28pt3_interp <- data.frame("depth" = depth_interp,
                                      "d13C_Ctot" = approx(d13C_BGveg_28pt3$depth,d13C_BGveg_28pt3$d13C_Ctot,depth_interp)$y)

d13C_BGveg_27pt3_interp <- data.frame("depth" = depth_interp,
                      "d13C_Ctot" = approx(d13C_BGveg_27pt3$depth,d13C_BGveg_27pt3$d13C_Ctot,depth_interp)$y)

d13C_rhizodeposits_29pt4_interp <- data.frame("depth" = depth_interp,
                      "d13C_Ctot" = approx(d13C_rhizodeposits_29pt4$depth,d13C_rhizodeposits_29pt4$d13C_Ctot,depth_interp)$y)

d13C_rhizodeposits_28pt4_interp <- data.frame("depth" = depth_interp,
                      "d13C_Ctot" = approx(d13C_rhizodeposits_28pt4$depth,d13C_rhizodeposits_28pt4$d13C_Ctot,depth_interp)$y)

f_CO2_uptake_0_interp <- data.frame("depth" = depth_interp,
                      "d13C_Ctot" = approx(f_CO2_uptake_0$depth,f_CO2_uptake_0$d13C_Ctot,depth_interp)$y)

f_CO2_uptake_0pt05_interp <- data.frame("depth" = depth_interp,
                      "d13C_Ctot" = approx(f_CO2_uptake_0pt05$depth,f_CO2_uptake_0pt05$d13C_Ctot,depth_interp)$y)

pCO2_fact_0pt0108_interp <- data.frame("depth" = depth_interp,
                      "d13C_Ctot" = approx(pCO2_fact_0pt0108$depth,pCO2_fact_0pt0108$d13C_Ctot,depth_interp)$y)

pCO2_fact_0pt0172_interp <- data.frame("depth" = depth_interp,
                      "d13C_Ctot" = approx(pCO2_fact_0pt0172$depth,pCO2_fact_0pt0172$d13C_Ctot,depth_interp)$y)

averageResult_interp <- data.frame("depth" = depth_interp,
                                       "d13C_Ctot" = approx(averageResult$depth,averageResult$d13C_Ctot,depth_interp)$y)

# --------------------------
# Global sensitivity analysis
# --------------------------

folder <- "SAFE output/20230801 PAWN ouput - isotopes"

PAWN_d13C_topsoil <- read.csv(paste(folder, "/", "Sensitivity_d13C_topsoil.csv", sep = ""))
PAWN_d13C_subsoil <- read.csv(paste(folder, "/", "Sensitivity_d13C_subsoil.csv", sep = ""))
PAWN_diff_d13C_topsoil_subsoil <- read.csv(paste(folder, "/", "Sensitivity_diff_d13C_topsoil_subsoil.csv", sep = ""))

# ------------------------------------------------------------------------------
# The plots are created - d13C
# ------------------------------------------------------------------------------

# -------------------------
# PAWN sensitivity analysis
# -------------------------

# The labels are added
labels <- c("d13C AGveg", "d13C BGveg", "d13C rhizodeposits", "fract CO2 uptake", "PCO2_fact", "dummy")
PAWN_d13C_topsoil$labels <- labels
PAWN_d13C_subsoil$labels <- labels
PAWN_diff_d13C_topsoil_subsoil$labels <- labels

# The dummy index is subtracted
PAWN_d13C_topsoil[1:5,1:3] <- PAWN_d13C_topsoil[1:5,1:3] - PAWN_d13C_topsoil[6,1]
PAWN_d13C_subsoil[1:5,1:3] <- PAWN_d13C_subsoil[1:5,1:3] - PAWN_d13C_subsoil[6,1]
PAWN_diff_d13C_topsoil_subsoil[1:5,1:3] <- PAWN_diff_d13C_topsoil_subsoil[1:5,1:3] - PAWN_diff_d13C_topsoil_subsoil[6,1]

# d13C - topsoil
p_topsoil <- PAWN_d13C_topsoil[1:5,] %>% 
      arrange(labels) %>%
      mutate(labels = factor(labels, levels = c("d13C AGveg", "d13C BGveg", "d13C rhizodeposits", "fract CO2 uptake", "PCO2_fact"))) %>%

      ggplot() +
      geom_pointrange(aes(x = labels, y = mid, ymin = lower, ymax = upper)) +
      scale_x_discrete(labels = c(expression(paste(delta^{13}, "C ") [italic(leaf)]), 
                                  expression(paste(delta^{13}, "C ") [italic(root)]), 
                                  expression(paste(delta^{13}, "C ") [italic(exudates)]), 
                                  expression(paste(alpha)), 
                                  expression(italic("S")))) +
      coord_flip() +
      scale_y_continuous(expand = c(0,0), limits = c(-0.05,1)) +
      labs(x = "", 
           y = "PAWN sensitivity index",
           title = expression(paste("Topsoil ", delta^{13}, "C"))) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text=element_text(size = 14))

annotate_figure(p_topsoil, fig.lab = "(A)")

# d13C - subsoil
p_subsoil <- PAWN_d13C_subsoil[1:5,] %>% 
  arrange(labels) %>%
  mutate(labels = factor(labels, levels = c("d13C AGveg", "d13C BGveg", "d13C rhizodeposits", "fract CO2 uptake", "PCO2_fact"))) %>%
  
  ggplot() +
  geom_pointrange(aes(x = labels, y = mid, ymin = lower, ymax = upper)) +
  scale_x_discrete(labels = c(expression(paste(delta^{13}, "C ") [italic(leaf)]), 
                              expression(paste(delta^{13}, "C ") [italic(root)]), 
                              expression(paste(delta^{13}, "C ") [italic(exudates)]), 
                              expression(paste(alpha)), 
                              expression(italic("S")))) +
  coord_flip() +
  scale_y_continuous(expand = c(0,0), limits = c(-0.02,0.55)) +
  labs(x = "", 
       y = "PAWN sensitivity index",
       title = expression(paste("Subsoil ", delta^{13}, "C"))) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size = 14))

# d13C - difference topsoil - subsoil
p_diff <- PAWN_diff_d13C_topsoil_subsoil[1:5,] %>% 
  arrange(labels) %>%
  mutate(labels = factor(labels, levels = c("d13C AGveg", "d13C BGveg", "d13C rhizodeposits", "fract CO2 uptake", "PCO2_fact"))) %>%
  
  ggplot() +
  geom_pointrange(aes(x = labels, y = mid, ymin = lower, ymax = upper)) +
  # geom_hline(data = PAWN_diff_d13C_topsoil_subsoil[6,], aes(yintercept = mid), linetype = "dashed") +
  scale_x_discrete(labels = c(expression(paste(delta^{13}, "C ") [italic(leaf)]), 
                              expression(paste(delta^{13}, "C ") [italic(root)]), 
                              expression(paste(delta^{13}, "C ") [italic(exudates)]), 
                              expression(paste(alpha)), 
                              expression(italic("S")))) +
  coord_flip() +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.55)) +
  labs(x = "", 
       y = "PAWN sensitivity index",
       title = expression(paste(Delta^{13}, "C topsoil - subsoil"))) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size = 14))

df_avg <- averageResult_interp
df_low <- d13C_AGveg_28pt9_interp
df_high <- d13C_AGveg_29pt9_interp
xMin <- -28
xMax <- -24
fontSize <- 14
titleName <- expression(paste(delta^{13}, "C ") [italic(leaf)])

createPlot_d13C <- function(df_avg, df_low, df_high, xMin, xMax, titleName, fontSize){
  
  polygonData <- data.frame("x" = c(df_low$d13C_Ctot[1:60], rev(df_high$d13C_Ctot[1:60])),
                            "y" = c(df_avg$depth[1:60], rev(df_avg$depth[1:60])))
  
  plotID <- ggplot() + 
    # The colored area between the extremes
    geom_polygon(data = polygonData, aes(x = x, y = y), alpha = .5, fill = "#3690c0") +
    # The average result
    geom_path(data = df_avg, aes(x = d13C_Ctot, y = depth)) +
    # The results for different input values
    geom_path(data = df_low, aes(x = d13C_Ctot, y = depth), linetype = "dashed") +
    geom_path(data = df_high, aes(x = d13C_Ctot, y = depth), linetype = "dashed") +
    
    scale_y_continuous(trans = "reverse", limits = c(0.6, 0), expand = c(0,0)) + 
    scale_x_continuous(position = "top", limits = c(xMin, xMax), expand = c(0,0)) + 
    labs(x=expression(paste(delta^{13}, "C (\u2030)"))) +
    annotate("text", x=-27.8, y=0.54, label= titleName, hjust = 0, size = 4) +
    annotate("segment", x = df_low$d13C_Ctot[5] - 0.05, xend = df_high$d13C_Ctot[5] + 0.05, y = df_low$depth[5], yend = df_low$depth[5], color = "black", arrow = arrow(ends = "both", length = (unit(0.2, "cm")))) +
    annotate("segment", x = df_low$d13C_Ctot[50] - 0.01, xend = df_high$d13C_Ctot[50] + 0.01, y = df_low$depth[50], yend = df_low$depth[50], color = "black", arrow = arrow(ends = "both", length = (unit(0.2, "cm")))) +
    annotate("text", x = df_low$d13C_Ctot[5] + 0.5, y = df_low$depth[5], label= paste("Delta == ", abs(round(as.numeric(df_low$d13C_Ctot[5] - df_high$d13C_Ctot[5]), digits = 2))), hjust = 0, size = 4, parse = TRUE) +
    annotate("text", x = df_low$d13C_Ctot[50] - 2.2, y = df_low$depth[50], label= paste("Delta == ", abs(round(as.numeric(df_low$d13C_Ctot[50] - df_high$d13C_Ctot[50]), digits = 2))), hjust = 0, size = 4, parse = TRUE) +
    
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text=element_text(size = fontSize),
          axis.title=element_text(size = fontSize),
          panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 1),
          axis.line = element_line(colour = "black", 
                                   size = 0),
          plot.margin = unit(c(0,0.5,0.5,0.5),"cm"))
  
  return(plotID)
  
}

# The plots are created

xMin <- -28
xMax <- -24
fontSize <- 14

yVal_text <- 0.575

df_low_d13C <- d13C_AGveg_28pt9_interp
df_high_d13C <- d13C_AGveg_29pt9_interp
p1 <- createPlot_d13C(averageResult_interp, df_low_d13C, df_high_d13C, xMin, xMax, expression(paste(delta^{13}, "C ") [italic(leaf)]), fontSize)
p1 <- p1 + labs(y = "Depth (m)") +
  annotate("text", x=-27.8, y=yVal_text, label= "[-29.9; -28.9]", hjust = 0, size = 3.5)
  
df_low_d13C <- d13C_BGveg_27pt3_interp
df_high_d13C <- d13C_BGveg_28pt3_interp
p2 <- createPlot_d13C(averageResult_interp, df_low_d13C, df_high_d13C, xMin, xMax, expression(paste(delta^{13}, "C ") [italic(root)]), fontSize)
p2 <- p2 + theme(axis.title.y = element_blank()) +
  annotate("text", x=-27.8, y=yVal_text, label= "[-28.3; -27.3]", hjust = 0, size = 3.5)

df_low_d13C <- d13C_rhizodeposits_28pt4_interp
df_high_d13C <- d13C_rhizodeposits_29pt4_interp
p3 <- createPlot_d13C(averageResult_interp, df_low_d13C, df_high_d13C, xMin, xMax, expression(paste(delta^{13}, "C ") [italic(exudates)]), fontSize)
p3 <- p3 + theme(axis.title.y = element_blank()) +
  annotate("text", x=-27.8, y=yVal_text, label= "[-29.4; -28.4]", hjust = 0, size = 3.5)

df_low_d13C <- f_CO2_uptake_0pt05_interp
df_high_d13C <- f_CO2_uptake_0_interp
p4 <- createPlot_d13C(averageResult_interp, df_low_d13C, df_high_d13C, xMin, xMax, expression(paste(alpha)), fontSize)
p4 <- p4 + theme(axis.title.y = element_blank()) +
  annotate("text", x=-27.8, y=yVal_text, label= "[0; 0.05]", hjust = 0, size = 3.5)

df_low_d13C <- pCO2_fact_0pt0172_interp
df_high_d13C <- pCO2_fact_0pt0108_interp
p5 <- createPlot_d13C(averageResult_interp, df_low_d13C, df_high_d13C, xMin, xMax, expression(italic("S")), fontSize)
p5 <- p5 + theme(axis.title.y = element_blank()) +
  annotate("text", x=-27.8, y=yVal_text, label= "[0.0108; 0.0172]", hjust = 0, size = 3.5)

# ------------------------------------------------------------------------------
# All plots for d13C are combined
# ------------------------------------------------------------------------------

lay <- rbind(c(rep(1,5), rep(2,5), rep(3,5)),
             c(4,4,4,5,5,5,6,6,6,7,7,7,8,8,8),
             c(4,4,4,5,5,5,6,6,6,7,7,7,8,8,8))

text.size <- 12

quartz(title = "Sensitivity analysis", width = 32*0.3937, height = 18*0.3937)
grid.arrange(annotate_figure(p_topsoil, fig.lab = "(A)", fig.lab.size = text.size), 
             annotate_figure(p_subsoil, fig.lab = "(B)", fig.lab.size = text.size),
             annotate_figure(p_diff, fig.lab = "(C)", fig.lab.size = text.size),
             annotate_figure(p5, fig.lab = "(D)", fig.lab.size = text.size),
             annotate_figure(p4, fig.lab = "(E)", fig.lab.size = text.size),
             annotate_figure(p3, fig.lab = "(F)", fig.lab.size = text.size),
             annotate_figure(p2, fig.lab = "(G)", fig.lab.size = text.size),
             annotate_figure(p1, fig.lab = "(H)", fig.lab.size = text.size),
             layout_matrix = lay)

# quartz.save("Figures for manuscript/Figure_13C_sensitivity.png", type="png")

# ggarrange(p_topsoil, p_subsoil, p_diff, p1, p2, p3, p4, p5, ncol = 5, labels = c("(A)", "(B)", "(C)", "(D)", "(E)"), layout_matrix = lay)

# The plot is saved
# g <- ggarrange(p1, p2, p3, p4, p5, ncol = 5, labels = c("(A)", "(B)", "(C)", "(D)", "(E)"))
# ggsave("Sensitivity_simple.png", g, device = "png", dpi = 300, units = "cm", width = 32, height = 12)

# ------------------------------------------------------------------------------
# The plots are created - d14C
# ------------------------------------------------------------------------------

# A function to create the plots
createPlot_d14C <- function(df_avg, df1, df2, xMin, xMax, titleName){
  
  polygonData <- data.frame("x" = c(df1$d14C_Ctot, rev(df2$d14C_Ctot)),
                            "y" = c(df_avg$depth, rev(df_avg$depth)))
  
  plotID <- ggplot() + 
    # The colored area between the extremes
    geom_polygon(data = polygonData, aes(x = x, y = y), alpha = .5, fill = "cadetblue2") +
    # The average result
    geom_path(data = df_avg, aes(x = d14C_Ctot, y = depth)) +
    # The results for different input values
    geom_path(data = df1, aes(x = d14C_Ctot, y = depth), linetype = "dashed") +
    geom_path(data = df2, aes(x = d14C_Ctot, y = depth), linetype = "dashed") +
    
    scale_y_continuous(trans = "reverse", limits = c(averageResult$depth[10], 0)) + 
    scale_x_continuous(position = "top", limits = c(xMin, xMax)) + 
    ggtitle(titleName) +
    labs(y="Depth (m)",
         x=expression(paste(Delta^{14}, "C (\u2030)"))) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(plotID)
  
}

# The plots are created

xMin <- -600
xMax <- 300

p1 <- createPlot_d14C(averageResult, d13C_AGveg_28pt9, d13C_AGveg_29pt9, xMin, xMax, "d13C AGveg")
p2 <- createPlot_d14C(averageResult, d13C_BGveg_28pt3, d13C_BGveg_27pt3, xMin, xMax, "d13C BGveg")
p3 <- createPlot_d14C(averageResult, diff_leaf_exudate_0, diff_leaf_exudate_1, xMin, xMax, "diff_leaf_exudate")
p4 <- createPlot_d14C(averageResult, f_CO2_uptake_0, f_CO2_uptake_0pt05, xMin, xMax, "f_CO2_uptake")
p5 <- createPlot_d14C(averageResult, pCO2_fact_0pt009, pCO2_fact_0pt019, xMin, xMax, "pCO2_fact")

quartz(title = "Effect of parameter values on isotope depth profiles", width = 8, height = 5)
grid.arrange(p1, p2, p3, p4, p5, ncol = 5)