# ----------------------------
# --- Weight Optimisation
# ============================
## Weight Estimations based on Raymer pp 460

## Aspect Ratio ======================================================================
# Effect of Aspect Ratio
# ARvarplot <- data.frame(AR = NA, MTOM = NA)
# for (AR in seq(18,100,2)){
#   W_dg_SI = ModifiedSecant(function(W_dg_SI) W_dg_SI - Weight_Estimate(WS, PW, W_dg_SI, composite = TRUE, iteration = TRUE),
#                            6000, 0.001, 0.01, positive = TRUE)
#   ARvarplot <- rbind(ARvarplot, data.frame(AR = AR, MTOM = W_dg_SI))
# }
# 
# ggplot(data = ARvarplot, aes(x = AR, y = MTOM)) + geom_point()

AR = 20


# Grid of data values
varWS <- seq(1500,3500, length.out = 11)
varPW <- seq(0, 30, length.out = 11)
weightoptim <- expand.grid(WS = varWS, PW = varPW)

weightoptim <- weightoptim %>%
  rowwise() %>%
  do(data.frame(
    # Previous WS and PW values
    .,
    # Determine MTOM
    MTOM = ModifiedSecant(
      function(W_dg_SI)
        W_dg_SI - Weight_Estimate(.$WS, .$PW, W_dg_SI, composite = TRUE, iteration = TRUE),
      6000, 0.001,0.01, positive = TRUE
    )
  )) %>%
  do(data.frame(
    ., 
    Weight_Estimate(.$WS, .$PW, .$MTOM, composite = TRUE, iteration = FALSE)[[1]]
  ))
weightoptim <- data.frame(weightoptim)

# Discrete Filled Contours
weightoptimdiscrete <- weightoptim
brks <- cut(weightoptim$MTOM, breaks = seq(3000, 20000, 500), dig.lab = 10)
brks <- gsub(",", " - ", brks, fixed = TRUE)
weightoptimdiscrete$MTOM <- gsub("\\(|\\]","",brks)

ggplot(data = weightoptimdiscrete, aes(x = WS, y = PW)) +
  geom_tile(aes(fill = MTOM))
# +
#   scale_fill_manual("MTOM",values=brewer.pal(14,"YlOrRd"))

MTOMlattice <- weightoptim %>%
  select(WS, PW, MTOM) %>%
  spread(PW, MTOM)
rownames(MTOMlattice) <- MTOMlattice$WS
MTOMlattice <- select(MTOMlattice, -WS)
asdf <- as.matrix(MTOMlattice)
filled.contour(varWS,varPW,asdf,nlevels=9,col=brewer.pal(9,"YlOrRd"))

# 
# points <- ggplot(data = weightoptim, aes(x = WS, y = PW)) +
#   geom_point(aes(colour = MTOM, size = 1/MTOM))

# MTOM contour plot
# Upside down parabolic shapes, centered around 2800
MTOMcontourplot <- ggplot(data = weightoptim, aes(x = WS, y = PW)) + 
  stat_contour(aes(z = MTOM, colour = ..level..)#breaks=c(c(seq(4000, 8000, 100), seq(7000, 20000, 500)))
               ) +
  geom_point(data = inp, aes(x = WS, y = P0/W, label = "Design")) +
  ggtitle("MTOM")
LabelledMTOMcontourplot <- direct.label(MTOMcontourplot, method = "top.pieces")

# Battery fraction contour plot
# Seems to be a mimium fraction around 2000, parabola shaping curves
BatteryFractioncontourplot <- ggplot(data = weightoptim, aes(x = WS, y = PW)) + 
  stat_contour(aes(z = batteries_Fraction, colour = ..level..)
    ) +
  ggtitle("Battery Fraction")
LabelledBatteryFractioncontourplot <- direct.label(BatteryFractioncontourplot, method = "top.pieces")

# Wing area (iterated solution)
# Decreases with increasing WS (duh) but more curvature at higher WS
Swingcontourplot <- ggplot(data = weightoptim, aes(x = WS, y = PW)) + 
  stat_contour(aes(z = S_wing, colour = ..level..),
               breaks = c(seq(5,23, 1), seq(25, 30, 2.5), seq(25,100,5))) +
  geom_point(data = inp, aes(x = WS, y = P0/W, label = "Design")) +
  ggtitle("Wing Area")
direct.label(Swingcontourplot, method = "top.pieces")

# Cd0 (iterated solution)
# Decreaes with WS since S_ref is increasing
Cd0contourplot <- ggplot(data = weightoptim, aes(x = WS, y = PW)) + 
  stat_contour(aes(z = Cd0, colour = ..level..)
  )
direct.label(Cd0contourplot, method = "top.pieces")

# Weight Savings (iterated solution)
# Decreaes with WS
Savingplot <- ggplot(data = weightoptim, aes(x = WS, y = PW))+ 
  stat_contour(aes(z = Saving_Total, colour = ..level..),
               breaks = c(seq(100, 150, 10), seq(150, 200, 25), seq(200, 1000, 50))
  ) +
  geom_point(data = inp, aes(x = WS, y = P0/W, label = "Design")) +
  ggtitle("Composite Weight Saving (Raymer Correlation)")
direct.label(Savingplot, method = "top.pieces")

# I should also plot out the tip chord and root chord!

ggplot(data = weightoptim) +
   geom_line(aes(x = S_wing, y = MTOM, group = PW, colour = PW))
ggplot(data = weightoptim) +
  geom_line(aes(x = S_wing, y = MTOM, group = WS, colour = WS))

ggplot(data = weightoptim) +
  geom_line(aes(x = WS, y = MTOM, colour = PW, group = PW))
ggplot(data = weightoptim) +
  geom_line(aes(x = PW, y = MTOM, colour = WS, group = WS))


carpetplot <- weightoptim %>%
  mutate(x = WS + MTOM/5)
# Carpet Plot (sort of)
ggplot(data = carpetplot, aes(x = x, y = MTOM, colour = MTOM)) +
  geom_path(aes(group = PW)) +
  geom_path(aes(group = WS)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())




LabelledMTOMcontourplot

## Constraint Plot ======================================================================
LabelledMTOMcontourplot + #ggplot(data = constraint, aes(x = WS, group = Clhls)) +
  geom_vline(data = constraint, aes(xintercept = WS_App, colour = "Landing")) +
  geom_line(data = constraint, aes(x = WS, group = Clhls, y = PW_TOFL, colour = "Takeoff")) +
  geom_line(data = constraint, aes(x = WS, group = Clhls, y = PW_Ceiling_Climb, colour = "Ceiling Climb")) +
  geom_line(data = constraint, aes(x = WS, group = Clhls, y = PW_Seg2_Climb, colour = "2nd Segment OEI")) +
  geom_line(data = constraint, aes(x = WS, group = Clhls, y = PW_Cruise_Climb, colour = "Cruise Climb")) +
  geom_vline(data = constraint, aes(xintercept = WS_Clstar, colour = "Clstar")) +
  geom_vline(data = constraint, aes(xintercept = WS_WbW0_Max, colour = "Wb/W0")) +
  geom_line(data = constraint, aes(x = WS, group = Clhls, y = PW_Cruise, colour = "Cruise")) +
  # Landing Labels
  geom_label(data = constraint, aes(x = WS_App, y = 12, label = sprintf("%0.2f", Clhls), 
                 colour = "Landing"), size = rel(3), show.legend = FALSE) + 
  geom_label(data = constraint, aes(x = WS_App, y = 12, label = sprintf("%0.2f", Clhls), 
                 colour = "Landing"), size = rel(3), show.legend = FALSE) + 
  # WbW0 Labels
  geom_label(data = constraint, aes(x = WS_WbW0_Max, y = 10, label = sprintf("%0.2f", WbW0_Max), 
                 colour = "Wb/W0"), size = rel(3), show.legend = FALSE) + 
  # Clstar Labels
  geom_label(data = constraint, aes(x = WS_Clstar, y = 14, label = sprintf("Cl*"), 
                 colour = "Clstar"), size = rel(3), show.legend = FALSE) + 
  # Cruise Labels
  geom_point(data = constraint, aes(x = WS, group = Clhls, y = PW_Cruise, colour = "Cruise")) +
  geom_label(data = constraint, aes(x = WS, group = Clhls, y = PW_Cruise, label = sprintf("%0.2f", Clcruise), 
                 colour = "Cruise"), size = rel(3), vjust = 1.5, show.legend = FALSE) +
  # Takeoff Labels
  geom_label(data = filter(constraint, WS == varWS[3]),
             aes(x = WS, y = PW_TOFL, label = sprintf("%0.2f", Clhls), 
                 colour = "Takeoff"), size = rel(3), show.legend = FALSE) + 
  # OEI Labels
  geom_label(data = filter(constraint, WS == varWS[2]),
             aes(x = WS, y = PW_Seg2_Climb, label = sprintf("%0.2f", Clhls), 
                 colour = "2nd Segment OEI"), size = rel(3), show.legend = FALSE) + 
  # Us
  geom_point(data = inp, aes(x = WS, y = P0/W, label = "Design")) +
  xlab("Wing Loading (N/m^2)") +
  ylab("Power Loading (W/N)") +
  ggtitle("Constraint Analysis") +
  geom_label(x = -Inf, y = Inf, vjust = 1.5, hjust = -0.1,
             label = "Constant Cd0, AR, e and eta")