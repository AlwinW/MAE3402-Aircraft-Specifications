# ----------------------------
# --- Weight Optimisation
# ============================
## Weight Estimations based on Raymer pp 460

## Constants ======================================================================
kg_to_lb = 2.20462
ft_to_m = 0.3048
Etatotal = inp$etamech*etaprop(inp$Vcruise)

Fudge_Factors = c(0.85, 0.83, 0.83, 0.90, 0.95, 0.95, 1, 1, 1, 1, 1, 1, 1, 1)

## Determine the Weight Estimate ======================================================================
Weight_Estimate <- function(WS, PW, W_dg_SI, composite = TRUE, iteration = TRUE,
                            AR = 20, BatteryFactor = 1.05, L_fusein = 12) {
  #--- Constants
  # Constant parameters for the aircraft
  # AR = 20
  AR_HT = 5.2
  AR_VT  = 1.76
  b_wing = 22.135/ft_to_m # ft, wing span
  e = 0.7993
  K_h = 0.11    # see pp461 includes hydraulics for flap
  lambda = 0.4
  L_Fuse = L_fusein/ft_to_m  # ft, fuselage length
  L_m = 0.7574/ft_to_m   # ft, length of main landing gear 
  L_n = 0.7574/ft_to_m   # ft, length of nose landing gear
  N_z = 3.5           # Ultimate load factor
  S_Eng = 2 * 2.281/ft_to_m^2  # ft^2 (both)
  S_Fuse = 55.551/ft_to_m^2    # ft^2
  t_c = 0.18            # Thickness to chord
  W_Fuse = 1.8/ft_to_m  # ft, Fuselage Width
  W_uav = 50*kg_to_lb   # lb, uninstalled avionics weight (800-1400 lb typ)
  # Constant aerodynamic properties
  L_D = 19 # (L/D)
  Mach = 0.25
  q = 0.5*17.56*10^-4*(0.25*1078)^2 # lb/ft^2
  q_SI = 1/2 * 0.9046291 * 82.09822^2 # Pa
  
  #--- Calculations
  # Directly Calculated from the Inputs
  S = W_dg_SI*9.8065 / WS        # m^2
  S_wing = S/(ft_to_m^2)         # ft^2
  P = PW * W_dg_SI * 9.8065      # kW
  W_dg = W_dg_SI*kg_to_lb        # lb
  # Engine sizing (engine + gearbox)
  W_en = P/(2 * 260e3) * (50 + 15) * kg_to_lb # lb, each engine
  # Tail sizing (FIX UP FOR VALUES TO MATCH DANIEL'S VALUES)
  L_t = 0.53*L_Fuse     # ft, wing to tail length
  S_HT = (0.07 * b_wing * S_wing)/L_t # ft^2 
  S_VT = (0.8 * 1.17456 / ft_to_m * S_wing)/L_t * 2 # ft^2
  # Cd0 estimate (FIX UP FOR WETTED VS THEORETICAL)
  S_wet = S_Fuse + S_Eng + 2 * S_wing + 2 * S_HT + 2 * S_VT - 10/(ft_to_m^2) # ft^2
  K = 1/(pi*AR*e)
  Cd0 = 0.0045 * S_wet / S_wing
  W_bat = (q_SI  * Cd0 / WS + K / q_SI * WS ) * BatteryFactor * 9.8065 / Etatotal * W_dg
  # Raymer Weights Correlation
  WeightsRow <- data.frame(
    wing = ((0.036*S_wing^0.758)*W_bat^0.0035*((AR)^0.6)*(q^0.006)*(lambda^0.04)*((100*(t_c))^-0.3)*((N_z*W_dg)^0.49)),
    htail = (0.016*(N_z*W_dg)^0.414)*(q^0.168)*(S_HT^0.896)*((100*(t_c))^-0.12)*(AR_HT)^0.043*(lambda^-0.02),
    vtail = 0.073*((N_z*W_dg)^0.376)*q^0.122*(S_VT^0.873)*((100*t_c)^-0.49)*((AR_VT)^0.357)*lambda^0.039,
    fuselage = 0.052*(S_Fuse^1.086)*((N_z*W_dg)^0.177)*(L_t^-0.051)*(L_D^-0.072)*(q^0.241),
    lgmain = (0.095*(N_z*W_dg)^0.768)*(L_m^0.409),
    lgnose = (0.125*(N_z*W_dg)^0.566)*(L_n^0.845),
    powerplant = 2.575*(W_en^0.922)*2,
    batteries = W_bat,
    payload = 6*120*kg_to_lb,
    controls = 0.053*L_Fuse^1.536*b_wing^0.371*(N_z*W_dg*10^-4)^0.8,
    hydraulics = K_h*(W_Fuse^0.8)*(Mach^0.5),
    avionics = 2.117*(W_uav)^0.933,
    ac_antiice = (0.265*W_dg^0.52)*(6^0.68)*((2.117*(W_uav)^0.933)^0.17)*(Mach^0.08),
    furnishings = 0.0582 * W_dg - 65
  )
  # Determine if advanced composites are used or not
  if (composite == TRUE)
    Fudge_Factors = c(0.85, 0.83, 0.83, 0.90, 0.95, 0.95, 1, 1, 1, 1, 1, 1, 1, 1)
  # OR ANOTHER CONDITION WITH PETE'S FACTORS INSTEAD?
  else
    Fudge_Factors = rep(1, 14)
  W_est = sum(WeightsRow*Fudge_Factors)/kg_to_lb
  # Output the result
  if (iteration == TRUE)
    # Quick one value output for convergence
    return(W_est)
  else{
    # Full output in a list
    # List of variables used
    Working_SI <- data.frame(
      S_wing = S,
      P = P,
      W_en = W_en / kg_to_lb,
      S_HT = S_HT * (ft_to_m)^2,
      S_VT = S_VT * (ft_to_m)^2,
      S_wet = S_wet * (ft_to_m)^2,
      Cd0 = Cd0,
      W_bat = W_bat / kg_to_lb
    )
    # Tabulated Weights - Imperial and SI
    Weights <- WeightsRow %>%
      gather(key = Part, value = General_Aviation) %>%
      data.frame(., Advanced_Composites = t(WeightsRow * Fudge_Factors))
    Total_Weights = data.frame(
      Part = "Total",
      General_Aviation = sum(Weights[, 2]), 
      Advanced_Composites = sum(Weights[, 3]))
    Weights <- rbind(Weights, Total_Weights)
    rownames(Weights) <- NULL
    Total_Weights_SI <- data.frame(
      Part = Total_Weights[, 1],
      Total_Weights[, 2:3] / kg_to_lb) %>%
      mutate(Saving = General_Aviation - Advanced_Composites)
    Weights_SI <- data.frame(
      Part = Weights[,1], 
      Weights[,2:3]/kg_to_lb)
    # Tabulated Weights Fraction
    Weights_Fraction <- data.frame(
      Part = Weights[,1], 
      General_Aviation = Weights[,2]/as.double(Total_Weights[2]), 
      Advanced_Composites = Weights[,3]/as.double(Total_Weights[3]))
    # Row data frame for cbind
    Rowsummary <- cbind(
      W_est = W_est,
      set_names(Total_Weights_SI[,2:4], paste(colnames(Total_Weights_SI[,2:4]),"Total", sep="_")),
      WeightsRow/kg_to_lb,
      set_names(WeightsRow/kg_to_lb/W_est, paste(colnames(WeightsRow),"Fraction", sep="_")),
      Working_SI
    )
    return(list(
      Rowsummary = Rowsummary,
      Weights = Weights_SI,
      Total_Weights_SI =  Total_Weights_SI,
      Weights_Fraction = Weights_Fraction,
      Working_SI = Working_SI
    ))
  }
}
## Determine the MTOW ======================================================================
WS = inp$WS
PW = inp$Pshaft/inp$W
W_dg_SI = ModifiedSecant(function(W_dg_SI) W_dg_SI - Weight_Estimate(WS, PW, W_dg_SI, composite = TRUE, iteration = TRUE), 
                         6000, 0.001, 0.01, positive = TRUE)
