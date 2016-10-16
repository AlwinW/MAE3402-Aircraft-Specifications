# ----------------------------
# --- Raymer Weight Correlations
# ============================

## Weight Estimations based on Raymer pp 460

# input: WS and PW and W

kg_to_lb = 2.20462
ft_to_m = 0.3048


S = W*9.8065 / WS
P = PW * W * 9.8065
W_en = P/(2 * 260e3) * 50 / 0.45 # SI (kg) to pound mass



S_wing = S/(ft_to_m^2)
W_in = W*kg_to_lb
Mach = 0.25
b_wing = 22.135/ft_to_m
lambda = 0.4
q = 0.5*17.56*10^-4*(0.25*1078)^2
AR = 20
t_c = 0.18
N_z = 3.5
L_t = 6.6/ft_to_m # wing to tail SI
L_m=0.8/ft_to_m # SI 
L_n = 0.8/ft_to_m # SI 
S_HT = 2.41/ft_to_m^2 # SI
AR_HT = 7.31
S_VT = 1.11/ft_to_m^2 # SI
AR_VT  = 1.77
S_Fuse = 51.5/ft_to_m^2 # SI
L_Fuse = 12/ft_to_m # SI Conversion
W_Fuse = 2.08/ft_to_m # Fuselage Width
L_D = 30 # (L/D)
K_h = 0.11 # see pp461 includes hydraulics for flap

WE = rep(0, 14)
WE[1] = ((0.036*S_wing^0.758)*((AR)^0.6)*(q^0.006)*(lambda^0.04)*((100*(t_c))^-0.3)*((N_z*W_in)^0.49))
WE[2] = (0.016*(N_z*W_in)^0.414)*(q^0.168)*(S_HT^0.896)*((100*(t_c))^-0.12)*(AR_HT)^0.043*(lambda^-0.02)
WE[3] = 0.073*((N_z*W_in)^0.376)*q^0.122*(S_VT^0.873)*((100*t_c)^-0.49)*((AR_VT)^0.357)*lambda^0.039
WE[4] = 0.052*(S_Fuse^1.086)*((N_z*W_in)^0.177)*(L_t^-0.051)*(L_D^-0.072)*(q^0.241)
WE[5] = (0.095*(N_z*W_in)^0.768)*(L_m^0.409)
WE[6] = (0.125*(N_z*W_in)^0.566)*(L_n^0.845)
WE[7] = 2.575*(W_en^0.922)*2 # Includes Prop + mounts
WE[8] = 2750*1.3*kg_to_lb
WE[9] = 6*120*kg_to_lb
WE[10] = 0.053*L_Fuse^1.536*b_wing^0.371*(N_z*W_in*10^-4)^0.8
WE[11] = K_h*(W_Fuse^0.8)*(Mach^0.5)
WE[12] = 50/0.45
WE[13] = (0.265*W_in^0.52)*(6^0.68)*((WE[12])^0.17)*(0.25^0.08)
WE[14] = 0.051*W_in - 65

W_out = sum(WE)/kg_to_lb
