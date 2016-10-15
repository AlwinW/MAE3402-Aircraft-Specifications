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
varWS <- seq(1500,3500, length.out = 21)
varPW <- seq(0, 30, length.out = 21)
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


# 
# points <- ggplot(data = weightoptim, aes(x = WS, y = PW)) +
#   geom_point(aes(colour = MTOM, size = 1/MTOM))

# MTOM contour plot
# Upside down parabolic shapes, centered around 2800
MTOMcontourplot <- ggplot(data = weightoptim, aes(x = WS, y = PW)) + 
  stat_contour(aes(z = MTOM, colour = ..level..)#breaks=c(c(seq(4000, 8000, 100), seq(7000, 20000, 500)))
               )
direct.label(MTOMcontourplot, method = "top.pieces")

# Battery fraction contour plot
# Seems to be a mimium fraction around 2000, parabola shaping curves
BatteryFractioncontourplot <- ggplot(data = weightoptim, aes(x = WS, y = PW)) + 
  stat_contour(aes(z = batteries_Fraction, colour = ..level..)
  )
direct.label(BatteryFractioncontourplot, method = "top.pieces")

# Wing area (iterated solution)
# Decreases with increasing WS (duh) but more curvature at higher WS
Swingcontourplot <- ggplot(data = weightoptim, aes(x = WS, y = PW)) + 
  stat_contour(aes(z = S_wing, colour = ..level..)
  )
direct.label(Swingcontourplot, method = "top.pieces")

Cd0contourplot <- ggplot(data = weightoptim, aes(x = WS, y = PW)) + 
  stat_contour(aes(z = Cd0, colour = ..level..)
  )
direct.label(Cd0contourplot, method = "top.pieces")

# Cd0 (iterated solution)
# Decreaes with WS since S_ref is increasing
Savingplot <- ggplot(data = weightoptim, aes(x = WS, y = PW)) + 
  stat_contour(aes(z = Saving_Total, colour = ..level..)
  )
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
