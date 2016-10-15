# ----------------------------
# --- Weight Optimisation
# ============================
## Weight Estimations based on Raymer pp 460

# Effect of Aspect Ratio
# ARvarplot <- data.frame(AR = NA, MTOM = NA)
# for (AR in seq(18,100,2)){
#   W_dg_SI = ModifiedSecant(function(W_dg_SI) W_dg_SI - Weight_Estimate(WS, PW, W_dg_SI, composite = TRUE, iteration = TRUE), 
#                            6000, 0.001, 0.01, positive = TRUE)
#   ARvarplot <- rbind(ARvarplot, data.frame(AR = AR, MTOM = W_dg_SI))
# }
# 
# ggplot(data = ARvarplot, aes(x = AR, y = MTOM)) + geom_point()


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

points <- ggplot(data = weightoptim, aes(x = WS, y = PW)) +
  geom_point(aes(colour = MTOM, size = 1/MTOM))

MTOMcontourplot <- ggplot(data = weightoptim, aes(x = WS, y = PW)) + 
  stat_contour(aes(z = MTOM, colour = ..level..)#, 
               # breaks=c(c(seq(4000, 8000, 100), seq(7000, 20000, 500)))
               )
direct.label(MTOMcontourplot, method = "top.pieces")

BatteryFractioncontourplot <- ggplot(data = weightoptim, aes(x = WS, y = PW)) + 
  stat_contour(aes(z = batteries_Fraction, colour = ..level..)
  )
direct.label(BatteryFractioncontourplot, method = "top.pieces")

Swingcontourplot <- ggplot(data = weightoptim, aes(x = WS, y = PW)) + 
  stat_contour(aes(z = S_wing, colour = ..level..)
  )
direct.label(Swingcontourplot, method = "top.pieces")
