# ----------------------------
# --- Weight Optimisation
# ============================
## Weight Estimations based on Raymer pp 460

ARvarplot <- data.frame(AR = NA, MTOM = NA)
for (AR in seq(18,100,2)){
  W_dg_SI = ModifiedSecant(function(W_dg_SI) W_dg_SI - Weight_Estimate(WS, PW, W_dg_SI, composite = TRUE, iteration = TRUE), 
                           6000, 0.001, 0.01, positive = TRUE)
  ARvarplot <- rbind(ARvarplot, data.frame(AR = AR, MTOM = W_dg_SI))
}

ggplot(data = ARvarplot, aes(x = AR, y = MTOM)) + geom_point()
