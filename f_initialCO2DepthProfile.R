model_initCO2DepthProfile <- function(Time, State, Pars, paramValues){
  
  with(as.list(c(State, Pars)), {
    
    # The state variables are stored in new variables
    CO2_12 <- State[startsWith(names(State), "CO2_12")]
    CO2_13 <- State[startsWith(names(State), "CO2_13")]
    CO2_14 <- State[startsWith(names(State), "CO2_14")]
  
    # Depth profiles of the production of 12CO2, 13CO2 and 14CO2 
    P_12C <- annualRootCO2Production_12C
    P_13C <- annualRootCO2Production_13C
    P_14C <- annualRootCO2Production_14C
    
    # Depth profiles of CO2 concentration of the atmosphere and soil
    conc_12C <- c(atmCO2_12C_MolPerCubicMeter, CO2_12, tail(CO2_12,1))
    conc_13C <- c(atmCO2_13C_MolPerCubicMeter, CO2_13, tail(CO2_13,1))
    conc_14C <- c(atmCO2_14C_MolPerCubicMeter, CO2_14, tail(CO2_14,1))
    
    # The depths are constructed
    depth <- c(0, paramValues@midDepth_layer, tail(paramValues@boundaryDepth,1))
    
    # dydt is initiated
    dydt = list()
    
    # --------------------------------------------------------------------------
    # The change in state variables is calculated following Goffin et al. (2014)
    # --------------------------------------------------------------------------
    
    # -----
    # 12CO2
    # -----
    
    #Harmonic mean of Ds
    Ds_full_12C <- c(D0_atm_year_12C, De_12C_year, tail(De_12C_year,1))
    endOf <- length(Ds_full_12C)
    Ds_top_12C <- harmonic.mean(rbind(Ds_full_12C[1:(endOf-2)], Ds_full_12C[2:(endOf-1)]))
    Ds_bot_12C <- harmonic.mean(rbind(Ds_full_12C[2:(endOf-1)], Ds_full_12C[3:endOf]))
    
    # The layer thickness
    endOf <- length(depth)
    dz_top <- depth[2:(endOf-1)] - depth[1:(endOf-2)]
    dz_bot <- depth[3:endOf] - depth[2:(endOf-1)]
    dz <- paramValues@layerThickness
    
    # dCO2 is defined
    dCO2_12C <- conc_12C
    endOf <- length(dCO2_12C)
    dCO2_top_12C <- dCO2_12C[1:(endOf-2)] - dCO2_12C[2:(endOf-1)]
    dCO2_bot_12C <- dCO2_12C[2:(endOf-1)] - dCO2_12C[3:endOf]
    
    # The fluxes through the top and bottom of every layer
    Ftop_12C <- -Ds_top_12C * (dCO2_top_12C / dz_top)
    Fbot_12C <- -Ds_bot_12C * (dCO2_bot_12C / dz_bot)
    
    # The change in CO2 concentration is calculated
    dydt$CO2_12 <- P_12C - ((Ftop_12C - Fbot_12C) / dz)
    
    # -----
    # 13CO2
    # -----
    
    #Harmonic mean of Ds
    Ds_full_13C <- c(D0_atm_year_13C, De_13C_year, tail(De_13C_year,1))
    endOf <- length(Ds_full_13C)
    Ds_top_13C <- harmonic.mean(rbind(Ds_full_13C[1:(endOf-2)], Ds_full_13C[2:(endOf-1)]))
    Ds_bot_13C <- harmonic.mean(rbind(Ds_full_13C[2:(endOf-1)], Ds_full_13C[3:endOf]))
    
    # The layer thickness
    endOf <- length(depth)
    dz_top <- depth[2:(endOf-1)] - depth[1:(endOf-2)]
    dz_bot <- depth[3:endOf] - depth[2:(endOf-1)]
    dz <- paramValues@layerThickness
    
    # dCO2 is defined
    dCO2_13C <- conc_13C
    endOf <- length(dCO2_13C)
    dCO2_top_13C <- dCO2_13C[1:(endOf-2)] - dCO2_13C[2:(endOf-1)]
    dCO2_bot_13C <- dCO2_13C[2:(endOf-1)] - dCO2_13C[3:endOf]
    
    # The fluxes through the top and bottom of every layer
    Ftop_13C <- -Ds_top_13C * (dCO2_top_13C / dz_top)
    Fbot_13C <- -Ds_bot_13C * (dCO2_bot_13C / dz_bot)
    
    # The change in CO2 concentration is calculated
    dydt$CO2_13 <- P_13C - ((Ftop_13C - Fbot_13C) / dz)
    
    # -----
    # 14CO2
    # -----
    
    #Harmonic mean of Ds
    Ds_full_14C <- c(D0_atm_year_14C, De_14C_year, tail(De_14C_year,1))
    endOf <- length(Ds_full_14C)
    Ds_top_14C <- harmonic.mean(rbind(Ds_full_14C[1:(endOf-2)], Ds_full_14C[2:(endOf-1)]))
    Ds_bot_14C <- harmonic.mean(rbind(Ds_full_14C[2:(endOf-1)], Ds_full_14C[3:endOf]))
    
    # The layer thickness
    endOf <- length(depth)
    dz_top <- depth[2:(endOf-1)] - depth[1:(endOf-2)]
    dz_bot <- depth[3:endOf] - depth[2:(endOf-1)]
    dz <- paramValues@layerThickness
    
    # dCO2 is defined
    dCO2_14C <- conc_14C
    endOf <- length(dCO2_14C)
    dCO2_top_14C <- dCO2_14C[1:(endOf-2)] - dCO2_14C[2:(endOf-1)]
    dCO2_bot_14C <- dCO2_14C[2:(endOf-1)] - dCO2_14C[3:endOf]
    
    # The fluxes through the top and bottom of every layer
    Ftop_14C <- -Ds_top_14C * (dCO2_top_14C / dz_top)
    Fbot_14C <- -Ds_bot_14C * (dCO2_bot_14C / dz_bot)
    
    # The change in CO2 concentration is calculated
    dydt$CO2_14 <- P_14C - ((Ftop_14C - Fbot_14C) / dz)
    
    dydt <- list(c(dydt$CO2_12, dydt$CO2_13, dydt$CO2_14))

    return(dydt)
  
  })
  
}