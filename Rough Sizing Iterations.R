#----------------------------
#--- Sizing Iterations
#============================

rootfindWbWo <- function(Mb, Sref, qinf, Cd0cor, K, Mnobatt, Etatotal, g_sl) {
  WbW0 = (qinf * Cd0cor/(((Mb + Mnobatt) * g_sl)/Sref) + K/qinf * (((Mb + Mnobatt) * g_sl)/Sref)) * g_sl* 1.1 / Etatotal - Mb/(Mb + Mnobatt)
}

# x = seq(3000, 6000, by = 100)
# plot(x, rootfindWbWo(x, 24.5, 3048.653, 0.0228938,  0.01924, 2725.182, 0.73, 9.8065))
# abline(h = 0, v = 4555.185, col = "lightgray", lty = 3)

ModifiedSecant(function(Mb) {
  rootfindWbWo(Mb, 24.5, 3048.653, 0.0228938,  0.01924, 2725.182, 0.73, 9.8065)},
  1000, 0.001, 0.001)


# Mass breakdown
mass_breakdown <- data.frame(
  wing = 566.1,
  htail = 18.8,
  vtail = 8.1,
  fuselage = 279,
  mainlg = 263,
  noselg = 59,
  engines = 180,
  batteries = 3575,
  payload = 720,
  controls = 122,
  hydraulics = 22.68,
  avionics = 50.4,
  airconn = 118,
  furnishing = 306
)

# Surface areas
Sref_breakdown <- data.frame(
    wing = 24.5,
    htail = 2.412,
    vtail = 1.11
  ) %>%
  mutate(
    hratio = htail/wing,
    vratio = vtail/wing,
    wingdens = mass_breakdown$wing/wing,
    htaildens = mass_breakdown$htail/htail,
    vtaildens = mass_breakdown$vtail/vtail)

# Mass summaries
mass_summary <- data.frame(
  MTOM = sum(mass_breakdown),
  fixed = sum(select(mass_breakdown, -wing, -htail, -vtail, -engines, -batteries, -payload)),
  payload = mass_breakdown$payload,
  wing = mass_breakdown$wing,
  tail = sum(select(mass_breakdown, htail, vtail)),
  batteries = mass_breakdown$batteries,
  engines = mass_breakdown$engines
)

# Vary the wing surface areas
varSrefwing <- seq(17, 45, by = 1)
# varSrefwing = 25
wingsurf <- data.frame(
    wingref = varSrefwing,
    Cfe = 0.0040, ## NOT 0.0045, adjusted to meed our 0.0191!!
    Swetfixed = 70) %>%
  mutate(
    htailref = wingref * Sref_breakdown$hratio,
    vtailref = wingref * Sref_breakdown$vratio,
    Swet = Swetfixed + (wingref + htailref + vtailref) * 2, # minus some interference
    Sref = wingref,
    SwetSref = Swet / Sref,
    Cd0cor = SwetSref * Cfe
    )

# Combine with inp parameters to find the weight and drag
wingopt <- cbind(inp, wingsurf) %>%
  select(-S) %>%
  mutate(h = AltCruise) %>%
  StandardAtomsphere(.) %>%
  mutate(
    Vinf = Mach * a,
    qinf = 1/2 * rho * Vinf^2,
    Mnobatt = 
      wingref * Sref_breakdown$wingdens + 
      htailref * Sref_breakdown$htaildens + 
      vtailref * Sref_breakdown$vtaildens + 
      mass_summary$fixed + 
      mass_summary$payload + 
      mass_summary$engines
  ) %>%
  rowwise() %>%
  do(data.frame(., Mbatt = ModifiedSecant(function(Mb) {
    rootfindWbWo(Mb, .$Sref, .$qinf, .$Cd0cor, .$K, .$Mnobatt, .$Etatotal, .$g_sl)},
    3000, 0.001, 0.001))) %>%
  ungroup() %>%
  mutate(WbW0 = (Mbatt)/(Mbatt + Mnobatt),
         MTOW = (Mbatt + Mnobatt)*g_sl,
         WS = (Mbatt + Mnobatt)/Sref * g_sl) %>%
  data.frame(.)

# Double check the solution
wingopt <- wingopt %>%
  mutate(
    WbW0_check = (qinf  * Cd0cor / WS + K / qinf * WS ) * 1.1 * g_sl / Etatotal,
    WbW0Cd0_check = (qinf  * Cd0cor / WS) * 1.1 * g_sl / Etatotal,
    WbW0K_check = (K / qinf * WS ) * 1.1 * g_sl / Etatotal,
    LHS_check = (qinf * Cd0cor/(WS) + K/qinf * (WS)) * g_sl* 1.1 / Etatotal,
    RHS_check = Mbatt/(Mbatt + Mnobatt)
  )




sampleplot <- select(wingopt, Sref, SwetSref, Cd0cor, MTOW, Swet, WS, WbW0) %>%
  gather(key = name, value = value, - Sref)

ggplot(data = sampleplot) +
  geom_line(aes(x = Sref, y = value)) +
  facet_wrap(~name, scales = "free_y")

# ggplot(data = wingopt) +
#   geom_line(aes(x = WS, y = WbW0))

# ggplot(data = wingopt, aes(x = WS)) +
#   geom_line(aes(y = WbW0_check)) + 
#   geom_line(aes(y = WbW0Cd0_check), linetype = 2) +
#   geom_line(aes(y = WbW0K_check), linetype = 2) +
#   ggtitle("Variable Cd0")
