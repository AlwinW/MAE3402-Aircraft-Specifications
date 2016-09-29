#----------------------------
#--- Constrain Analysis
#============================

varWS = seq(1000,3000, by = 250)
varClhls = seq(0.9,1.2, by = 0.3)

constraint <- inp %>%
  select(-S, -b, -m, -W, -P0eng, -P0) %>%
  mutate(h = AltCruise) %>%
  StandardAtomsphere(.)

constraint <-  RepeatRows(constraint, length(varClhls))
constraint$Clhls <- varClhls
constraint$Clflaps <- varClhls

constraint <-  RepeatRows(constraint, length(varWS))
constraint$WS <- varWS

constraint <- constraint %>%
  mutate(
    ClTO = Clclean + Clflaps,
    ClLD = Clclean + Clhls,
    Vcruise = Mach * a,
    qinf = 1/2 * rho * Vcruise^2,
## Landinng Approach Speed ======================================================================
    WS_App = (ClTO) * (1/2 * rho_sl) * (Vappmax / 1.3) ^ 2,
## Takeoff Field Length ======================================================================
    TOa = WS^2 * (0.225/ClTO^2),
    TOb = (11.8/ClTO) * WS,
    TOc = - Srun,
    PW_TOFL = 2*TOa/(-TOb + sqrt(TOb^2 - 4*TOa*TOc)),
## Cruise Speed ======================================================================
    Clcruise = WS / qinf,
    PW_Cruise = (Vcruise/Etaprop) * (Cd0/Clcruise + K*Clcruise),
## Climb at Ceilinng ======================================================================
    rho_ceil = 0.8491284, # DON'T HARD CODE IT LATER
    PW_Ceiling_Climb = ClimbCeil/Etaprop + (2/(Etaprop * rho_ceil)) * sqrt((K * WS)/(3*Cd0)) * (1.155 * sqrt(4*Cd0*K)),
## Climb OEI 2nd Segment ======================================================================
    c = sqrt(1.2^2 * 2 / (rho_sl * ClTO)),
    Clseg2 = ClTO / 1.2^2,
    Cdseg2 = Cd0G + K * Clseg2^2,
    PW_Seg2_Climb = (2 * (PerGrad2Seg/100 + Cdseg2/Clseg2) / (Etaprop/c)) * sqrt(WS),
## Climb at Cruise ======================================================================
    PW_Cruise_Climb = ClimbCruise/Etaprop + (2/(Etaprop * rho)) * sqrt((K * WS)/(3*Cd0)) * (1.155 * sqrt(4*Cd0*K)),
## Fly Near Clstar ======================================================================
    WS_Clstar = qinf * sqrt(Cd0/K),
## Empty Weight Fraction ======================================================================
    WbWe = (qinf  * Cd0 / WS + K / qinf * WS ) * 1.1 * g_sl / Etatotal,
    WbWe_Max = 0.6,
    WS_WbWe_Max = (WbWe_Max * qinf - sqrt(qinf^2 * (WbWe_Max^2 - 4*Cd0*(1.1*g_sl/Etatotal)^2*K))) / (2*(1.1*g_sl/Etatotal)*K)
  )

print(constraint)

ggplot(data = constraint, aes(x = WS)) +
  geom_vline(aes(xintercept = WS_App)) +
  geom_vline(aes(xintercept = WS_Clstar)) +
  geom_vline(aes(xintercept = WS_WbWe_Max)) +
  geom_line(aes(y = PW_TOFL)) +
  geom_line(aes(y = PW_Cruise)) +
  geom_line(aes(y = PW_Ceiling_Climb)) +
  geom_line(aes(y = PW_Seg2_Climb)) +
  geom_line(aes(y = PW_Cruise_Climb)) + 
  facet_grid(~Clhls)

# ggplot(data = constraint, aes(x = WS)) +
#   geom_line(aes(y = WbWe)) +
#   facet_grid(~Clhls)