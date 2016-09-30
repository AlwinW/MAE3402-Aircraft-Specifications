#----------------------------
#--- Sizing Iterations
#============================

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
varSrefwing <- seq(5, 50, by = 1)
# varSrefwing = 24.5
wingsurf <- data.frame(
    wingref = varSrefwing,
    Cfe = 0.0045,
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
  do(data.frame(., Mbatt = ModifiedSecant(function(Wb) {
    (.$qinf * .$Cd0cor/((Wb + .$Mnobatt)/.$Sref) + .$K/.$qinf * ((Wb + .$Mnobatt)/.$Sref)) * 1.1 / .$Etatotal - Wb/(Wb + .$Mnobatt)},
    3000, 0.001, 0.1))) %>%
  ungroup() %>%
  mutate(WbW0 = (Mbatt)/(Mbatt + Mnobatt),
         WS = (Mbatt + Mnobatt)/Sref * g_sl) %>%
  data.frame(.)

sampleplot <- select(wingopt, Sref, SwetSref, Cd0cor, Mnobatt, Mbatt, WbW0, WS) %>%
  gather(key = name, value = value, - Sref)

ggplot(data = sampleplot) +
  geom_line(aes(x = Sref, y = value)) +
  facet_wrap(~name, scales = "free_y")


