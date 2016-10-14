# ----------------------------
# --- Raymer Weight Correlations
# ============================

## Weight Estimations based on Raymer pp 460

## Constants ======================================================================
kg_to_lb = 2.20462
ft_to_m = 0.3048
Etatotal = 0.80
BatteryFactor = 1.0

# Constant parameters for the aircraft
AR = 20
b_wing = 22.135/ft_to_m # ft, wing span
K_h = 0.11    # see pp461 includes hydraulics for flap
lambda = 0.4
K = 0.01991
L_Fuse = 12/ft_to_m # ft, fuselage length
L_t = 6/ft_to_m     # ft, wing to tail length
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

# Test wingloading
WS = 2507.8
PW = 20

## Iterate out the weight ======================================================================
# Initial weight guess
W_dg = 8000*kg_to_lb


# Directly Calculated from the Inputs
S = W_dg*9.8065 / WS       # m^2
S_wing = S/(ft_to_m^2)  # ft^2
P = PW * W_dg * 9.8065     # kW

# Engine sizing
W_en = P/(2 * 260e3) * 50 * kg_to_lb # lb, each engine

# Tail sizing (FIX UP FOR VALUES TO MATCH DANIEL'S VALUES)
S_HT = (0.07 * b_wing * S_wing)/L_t  # ft^2 
AR_HT = 5.2
S_VT = (0.8 * 0.95 / ft_to_m * S_wing)/L_t * 2  # ft^2
AR_VT  = 1.76

# Cd0 estimate (FIX UP FOR WETTED VS THEORETICAL)
S_wet = S_Fuse + S_Eng + 2 * S_wing + 2 * S_HT + 2 * S_VT
Cd0 = 0.0045 * S_wet / S_wing
W_bat = (q_SI  * Cd0 / WS + K / q_SI * WS ) * BatteryFactor * 9.8065 / Etatotal * W_dg

Weights <- data.frame(
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

# WE = rep(0, 14)
# WE[1] = ((0.036*S_wing^0.758)*((AR)^0.6)*(q^0.006)*(lambda^0.04)*((100*(t_c))^-0.3)*((N_z*W_dg)^0.49))
# WE[2] = (0.016*(N_z*W_dg)^0.414)*(q^0.168)*(S_HT^0.896)*((100*(t_c))^-0.12)*(AR_HT)^0.043*(lambda^-0.02)
# WE[3] = 0.073*((N_z*W_dg)^0.376)*q^0.122*(S_VT^0.873)*((100*t_c)^-0.49)*((AR_VT)^0.357)*lambda^0.039
# WE[4] = 0.052*(S_Fuse^1.086)*((N_z*W_dg)^0.177)*(L_t^-0.051)*(L_D^-0.072)*(q^0.241)
# WE[5] = (0.095*(N_z*W_dg)^0.768)*(L_m^0.409)
# WE[6] = (0.125*(N_z*W_dg)^0.566)*(L_n^0.845)
# WE[7] = 2.575*(W_en^0.922)*2 # Includes Prop + mounts
# WE[8] = 2750*1.3*kg_to_lb
# WE[9] = 6*120*kg_to_lb
# WE[10] = 0.053*L_Fuse^1.536*b_wing^0.371*(N_z*W_dg*10^-4)^0.8
# WE[11] = K_h*(W_Fuse^0.8)*(Mach^0.5)
# WE[12] = 50/0.45
# WE[13] = (0.265*W_dg^0.52)*(6^0.68)*((WE[12])^0.17)*(0.25^0.08)
# WE[14] = 0.051*W_dg - 65

W_out = sum(Weights)/kg_to_lb
