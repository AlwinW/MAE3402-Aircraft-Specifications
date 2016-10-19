# ----------------------------
# --- Weight Optimisation Graphs
# ============================

## Contour Plots ======================================================================
# Grid of data values
varWS <- seq(1500,5500, length.out = 21)
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


MTOMlattice <- weightoptim %>% 
  mutate(MTOM = ifelse(MTOM > 12000, 12000, MTOM)) %>%
  select(WS, PW, MTOM) %>% 
  spread(PW, MTOM) 
rownames(MTOMlattice) <- MTOMlattice$WS 
MTOMlattice <- select(MTOMlattice, -WS) 
asdf <- as.matrix(MTOMlattice) 

filled.contour(varWS, varPW, asdf, nlevels = 9, col=brewer.pal(9,"YlOrRd"),
               xlab="Wing Loading (WS)",ylab="Power Loading (PW)",
               main = "Maximum Take Off Mass (MTOM)",
               # plot.title = "MTOM",
               plot.axes = {
                 axis(1); axis(2);
                 contour(varWS, varPW, asdf, 
                     drawlabels = TRUE, add = TRUE);
                 points(2000, 10);
                 lines(x = c(2000, 3000), y = c(10, 20))
                 })

filled.contour(varWS, varPW, asdf, nlevels=9, col=brewer.pal(9,"YlOrRd"),
               plot.axes = {axis(1); axis(2);
                 contour(x=varWS,y=varPW,z=asdf,nlevels=9,labcex=1.5,
                         col=pt.color,lwd=1.5,add=TRUE, labels="      ", method="flattest"
                 );
                 contour(x=varWS,y=varPW,z=asdf,nlevels=9,lwd=1.5,labcex=1.5,
                         lty=0,col=pt.color,add=TRUE, method="flattest"
                 );
                 points(2000, 10);
                 }) 


## 3D Plots ======================================================================
varWS <- seq(1500,3500, length.out = 11)
varPW <- seq(0, 30, length.out = 11)
varAR <- seq(15, 45, length.out = 7)
weightoptim3D <- expand.grid(WS = varWS, PW = varPW)
weightoptim3D <- RepeatRows(weightoptim3D, varAR)
weightoptim3D$AR <- varAR

weightoptim3D <- weightoptim3D %>%
  rowwise() %>%
  do(data.frame(
    # Previous WS and PW values
    .,
    # Determine MTOM
    MTOM = ModifiedSecant(
      function(W_dg_SI)
        W_dg_SI - Weight_Estimate(.$WS, .$PW, W_dg_SI, AR =.$AR, composite = TRUE, iteration = TRUE),
      6000, 0.001,0.01, positive = TRUE
    )
  ))
weightoptim3D <- data.frame(weightoptim3D)

plot_ly(data = weightoptim3D, x = ~WS, y = ~PW, z = ~AR, color = ~MTOM, alpha = 0.3)


