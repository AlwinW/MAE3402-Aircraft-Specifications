#----------------------------
#--- Constrain Analysis
#============================

# NOTE:
# I probably should have set it up like this:
#   constraint PW curves
#   constraint WS curves

varWS = seq(1500,3500, by = 250)
varClhls = seq(0.9,1.3, by = 0.1)
varWbW0_Max = seq(0.6,0.64, by = 0.01)

constraint <- inp %>%
  select(-S, -b, -m, -W, -P0eng, -P0) %>%
  mutate(h = AltCruise,
         Etaprop = 0.77,
         Etatotal = 0.73,
         BatteryFactor = 1.0) %>%
  StandardAtomsphere(.)

constraint <-  RepeatRows(constraint, length(varClhls))
constraint$Clhls <- varClhls
constraint$Clflaps <- varClhls

constraint <-  RepeatRows(constraint, length(varWS))
constraint$WS <- varWS

constraint <-  RepeatRows(constraint, length(varWbW0_Max))
constraint$WbW0_Max <- varWbW0_Max

constraint <- constraint %>%
  mutate(
    ClTO = Clclean + Clflaps,
    ClLD = Clclean + Clhls,
    Vcruise = Mach * a,
    qinf = 1/2 * rho * Vcruise^2,
## Landing Approach Speed ======================================================================
    WS_App = (ClTO) * (1/2 * rho_sl) * (Vappmax / 1.3) ^ 2,
## Takeoff Field Length ======================================================================
    TOa = WS^2 * (0.255/ClTO^2),
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
    PW_Seg2_Climb = 2 *((PerGrad2Seg/100 + Cdseg2/Clseg2)  * sqrt(WS) / (Etaprop/c)),
## Climb at Cruise ======================================================================
    PW_Cruise_Climb = ClimbCruise/Etaprop + (2/(Etaprop * rho)) * sqrt((K * WS)/(3*Cd0)) * (1.155 * sqrt(4*Cd0*K)),
## Fly Near Clstar ======================================================================
    Clstar = sqrt(Cd0/K),
    WS_Clstar = qinf * Clstar,
## Empty Weight Fraction ======================================================================
    WbW0 = (qinf  * Cd0 / WS + K / qinf * WS ) * BatteryFactor * g_sl / Etatotal,
    WbW0Cd0 = (qinf  * Cd0 / WS) * BatteryFactor * g_sl / Etatotal,
    WbW0K = (K / qinf * WS ) * BatteryFactor * g_sl / Etatotal,
    WS_WbW0_Max = (WbW0_Max * qinf - sqrt(qinf^2 * (WbW0_Max^2 - 4*Cd0*(BatteryFactor*g_sl/Etatotal)^2*K))) / (2*(BatteryFactor*g_sl/Etatotal)*K)
  )

# print(constraint)





## Plot ======================================================================
ggplot(data = constraint, aes(x = WS, group = Clhls)) +
  geom_vline(aes(xintercept = WS_App, colour = "Landing")) +
  geom_line(aes(y = PW_TOFL, colour = "Takeoff")) +
  geom_line(aes(y = PW_Ceiling_Climb, colour = "Ceiling Climb")) +
  geom_line(aes(y = PW_Seg2_Climb, colour = "2nd Segment OEI")) +
  geom_line(aes(y = PW_Cruise_Climb, colour = "Cruise Climb")) +
  geom_vline(aes(xintercept = WS_Clstar, colour = "Clstar")) +
  geom_vline(aes(xintercept = WS_WbW0_Max, colour = "Wb/W0")) +
  geom_line(aes(y = PW_Cruise, colour = "Cruise")) +
  # Landing Labels
  geom_label(aes(x = WS_App, y = 12, label = sprintf("%0.2f", Clhls), 
                 colour = "Landing"), size = rel(3), show.legend = FALSE) + 
  geom_label(aes(x = WS_App, y = 12, label = sprintf("%0.2f", Clhls), 
                 colour = "Landing"), size = rel(3), show.legend = FALSE) + 
  # WbW0 Labels
  geom_label(aes(x = WS_WbW0_Max, y = 11, label = sprintf("%0.2f", WbW0_Max), 
                 colour = "Wb/W0"), size = rel(3), show.legend = FALSE) + 
  # Clstar Labels
  geom_label(aes(x = WS_Clstar, y = 13, label = sprintf("Cl*"), 
                 colour = "Clstar"), size = rel(3), show.legend = FALSE) + 
  # Cruise Labels
  geom_point(aes(y = PW_Cruise, colour = "Cruise")) +
  geom_label(aes(y = PW_Cruise, label = sprintf("%0.2f", Clcruise), 
                 colour = "Cruise"), size = rel(3), vjust = 1.5, show.legend = FALSE) +
  # Takeoff Labels
  geom_label(data = filter(constraint, WS == varWS[2]),
             aes(x = WS, y = PW_TOFL, label = sprintf("%0.2f", Clhls), 
               colour = "Takeoff"), size = rel(3), show.legend = FALSE) + 
  # Us
  geom_point(data = inp, aes(x = WS, y = P0/W, label = "Design")) +
  xlab("Wing Loading") +
  ylab("")

# ggplot(data = constraint, aes(x = WS, group = Clhls)) +
#   geom_line(aes(y = WbW0)) + 
#   geom_line(aes(y = WbW0Cd0), linetype = 2) +
#   geom_line(aes(y = WbW0K), linetype = 2) +
#   ggtitle("Constant Cd0, S values")
