##Weight Estimations based on Raymer pp 460

S_wing = 245/(03048^2)
Mach = 025
b_wing = 22135/03048 #SI conversion to feet
lambda = 04
q = 05*1756*10^-4*(025*1078)^2
AR = 20
t_c = 018
N_z = 35
W = 1450642 #pound mass
L_t = 66/03048 #wing to tail SI
L_m=08/03048 #SI 
L_n = 08/03048 #SI 
S_HT = 241/03048^2 #SI
AR_HT = 731
S_VT = 111/03048^2 #SI
AR_VT  = 177
S_Fuse = 515/03048^2 #SI
L_Fuse = 12/03048 #SI Conversion
W_Fuse = 208/03048 #Fuselage Width
L_D = 30 #(L/D)
W_en = 50/045 #SI (kg) to pound mass
K_h = 011 #see pp461 includes hydraulics for flaps

WE = rep(0, 14)
WE[1] = ((0.036*S_wing^0.758)*((AR)^0.6)*(q^0.006)*(lambda^0.04)*((100*(t_c))^-0.3)*((N_z*W)^0.49))

WE[2] = (0.016*(N_z*W)^0.414)*(q^0.168)*(S_HT^0.896)*((100*(t_c))^-0.12)*(AR_HT)^0.043*(lambda^-0.02)

WE[3] = 0.073*((N_z*W)^0.376)*q^0.122*(S_VT^0.873)*((100*t_c)^-0.49)*((AR_VT)^0.357)*lambda^0.039

WE[4] = 0.052*(S_Fuse^1.086)*((N_z*W)^0.177)*(L_t^-0.051)*(L_D^-0.072)*(q^0.241)

WE[5] = (0.095*(N_z*W)^0.768)*(L_m^0.409)

WE[6] = (0.125*(N_z*W)^0.566)*(L_n^0.845)

WE[7] = 2.575*(W_en^0.922)*2 #Includes Prop + mounts

WE[8] = 2750*1.3*2.20462

WE[9] = 6*120*2.20462

WE[10] = 0.053*L_Fuse^1.536*b_wing^0.371*(N_z*W*10^-4)^0.8

WE[11] = K_h*(W_Fuse^0.8)*(Mach^0.5)

WE[12] = 50/0.45

WE[13] = (0.265*W^0.52)*(6^0.68)*((WE[12])^0.17)*(0.25^0.08)

WE[14] = 0.051*W - 65


# matrix = {'wing' 'HT' 'VT' 'Fuse' 'LGM' 'LGN' 'Engine' 'Battery', 'Payload','Flight Controls','Hydraulics','Avionics','AirCon + anti-ice','furnishing'}
Weight_Tot = sum(WE)
Weight_Tot_SI = Weight_Tot * 0.453592

sprintf('Total mass (kg) is %f',Weight_Tot_SI)
