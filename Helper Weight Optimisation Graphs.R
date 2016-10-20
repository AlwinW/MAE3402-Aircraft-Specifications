# ----------------------------
# --- Weight Optimisation Graphs
# ============================

## Contour Plots ======================================================================
# Grid of data values
varWS <- seq(1500,3500, length.out = 21)
varPW <- seq(10, 30, length.out = 21)
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

filled.contour(varWS, varPW, asdf, nlevels = 8, col=brewer.pal(9,"YlOrRd"),
               xlab="Wing Loading (WS)",ylab="Power Loading (PW)",
               main = "Maximum Take Off Mass (MTOM)",
               # plot.title = "MTOM",
               plot.axes = {
                 axis(1); axis(2);
                 contour(varWS, varPW, asdf, 
                     drawlabels = TRUE, add = TRUE);
                 points(2000, 10);
                 polygon(x = c(xaxis$WS, 3500, 3500, 0), y = c(xaxis$PW, 30, 0, 0), 
                         density = 3, col = "grey40");
                 lines(x = xaxis$WS, y = xaxis$PW, lwd = 2);
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

## Projected Plots ======================================================================
PW_Seg2_Climb <- function(inp, WSinp, PWval = FALSE, all = FALSE) {
  constraint <- RepeatRows(inp, length(WSinp))
  constraint$WS <- WSinp
  constraint <- constraint %>%
    mutate(
      c = sqrt(1.2^2 * 2 / (rho_sl * ClTO)),
      Clseg2 = ClTO / 1.2^2,
      Cdseg2 = Cd0G + K * Clseg2^2,
      PW_Seg2_Climb = 2 *((PerGrad2Seg/100 + Cdseg2/Clseg2)  * sqrt(WS) / (EtapropG/c))
    )
  if (PWval == TRUE)
    return(constraint$PW_Seg2_Climb)
  else if (all == TRUE)
    return(contraint)
}

PW_TOFL <- function(inp, WSinp, PWval = FALSE, all = FALSE) {
  constraint <- RepeatRows(inp, length(WSinp))
  constraint$WS <- WSinp
  constraint <- constraint %>%
    mutate(
      TOa = WS^2 * (0.255/ClTO^2),
      TOb = (11.8/ClTO) * WS,
      TOc = - Srun,
      PW_TOFL = 2*TOa/(-TOb + sqrt(TOb^2 - 4*TOa*TOc))   * 0.95 # INTEGRAL METHOD BETTER
    )
  if (PWval == TRUE)
    return(constraint$PW_TOFL)
  else if (all == TRUE)
    return(contraint)
}

WS_App <- function(inp, WSval = FALSE, all = FALSE) {
  constraint <- inp %>%
    mutate(
      WS_App = (ClTO) * (1/2 * rho_sl) * (Vappmax / 1.3) ^ 2
    )
  if (WSval == TRUE)
    return(constraint$WS_App)
  else if (all == TRUE)
    return(contraint)
}

ProjectedxAxis <- function(inp, WSlower, WSupper, PWupper) {
  intersect <- inp %>%
    select(-S, -b, -m, -W) %>%
    mutate(h = AltCruise,
           Etaprop = etaprop(inp$Vcruise),
           Etatotal = inp$etamech*Etaprop,
           BatteryFactor = 1.03,
           Cd0G = inp$Cd0clean + inp$Cd0lg + inp$Cd0flaps + inp$Cdiflaps + inp$Cd0propfea, # feathered, 0.205 unfeathered
           EtapropG = etaprop(1.2 * inp$VsTO)) %>%
    StandardAtomsphere(.) %>%
    mutate(
      ClTO = Clclean + Clflaps,
      ClLD = Clclean + Clhls,
      qinf = 1/2 * rho * Vcruise^2)
  
  Seg2_TOFL_WS <- ModifiedSecant(
    function(WS) PW_Seg2_Climb(intersect, WS, PWval = TRUE) - PW_TOFL(intersect, WS, PWval = TRUE),
    xr = 2000, del = 0.001, toler = 0.01, positive = TRUE
  )
  
  App_WS = WS_App(intersect, WSval = TRUE)
  Seg2_App_PW = PW_Seg2_Climb(intersect, WS = App_WS, PWval = TRUE)
  TOFL_App_PW = PW_TOFL(intersect, WS = App_WS, PWval = TRUE)
  ## DON'T WORRY ABOUT GENERAL CASES ANY MORE, JUST CODE AS IF Seg2_TOFL_WS < App_WS
  varWS = seq(WSlower, WSupper, by = 50)
  
  xaxis1 <- c(varWS[varWS < Seg2_TOFL_WS], Seg2_TOFL_WS)
  xaxis1 <- data.frame(
    WS = xaxis1,
    PW = PW_Seg2_Climb(intersect, WSinp = xaxis1, PWval = TRUE)
  )
  xaxis2 <- c(varWS[varWS > Seg2_TOFL_WS & varWS < App_WS], App_WS)
  xaxis2 <- data.frame(
    WS = xaxis2,
    PW = PW_TOFL(intersect, WSinp = xaxis2, PWval = TRUE)
  )
  varPW = seq(0, PWupper, 0.5)
  
  xaxis3 <- c(PW_TOFL(intersect, WS = App_WS, PWval = TRUE), varPW[varPW > PW_TOFL(intersect, WS = App_WS, PWval = TRUE)])
  xaxis3 <- data.frame(
    WS = App_WS,
    PW = xaxis3
  )
  
  xaxis <- cbind(
    pseudox = seq(1:nrow(xaxis)),
    rbind(xaxis1, xaxis2, xaxis3)
    )
  return(xaxis)
}

ProjectedGraphs <- function(xaxis, inp, AR = 20, L_fusein = 12) {
  projectedgraphs <- RepeatRows(inp, xaxis)
  projectedgraphs$x <- xaxis$x
  projectedgraphs$WS <- xaxis$WS
  projectedgraphs$PW <- xaxis$PW
  
  projectedgraphs <- RepeatRows(projectedgraphs, length(AR))
  projectedgraphs$AR <- AR
  
  projectedgraphs <- RepeatRows(projectedgraphs, length(L_fusein))
  projectedgraphs$L_fusein <- L_fusein
  
  projectedgraphs <- projectedgraphs %>%
    rowwise() %>%
    do(data.frame(
      # Previous WS and PW values
      .,
      # Determine MTOM
      MTOM = ModifiedSecant(
        function(W_dg_SI)
          W_dg_SI - Weight_Estimate(.$WS, .$PW, W_dg_SI, AR = .$AR, L_fusein = .$L_fusein,
                                    composite = TRUE, iteration = TRUE),
        6000, 0.001,0.01, positive = TRUE
      )
    )) %>%
    do(data.frame(
      ., 
      Weight_Estimate(.$WS, .$PW, .$MTOM, composite = TRUE, iteration = FALSE)[[1]]
    ))
  projectedgraphs <- data.frame(projectedgraphs)
  return(projectedgraphs)
}

## Effect of AR ======================================================================
varAR = seq(15, 82, 1)

projectedgraphs <- ProjectedGraphs(xaxis, inp, AR = varAR)
projectedgraphs <- mutate(projectedgraphs, AR = factor(AR, levels = varAR, ordered = TRUE))

labels <- select(projectedgraphs, x, WS, PW) %>% distinct(.keep_all = TRUE) %>% arrange(x) %>%
  mutate(label = paste(round(WS,1), round(PW,1), sep="\n"))
labels <- rbind(data.frame(x = 0, WS = NA, PW = NA, label = "WS:\nPW:"), labels)

ggplot(projectedgraphs, aes(x = x, y = MTOM, group = AR, colour = AR)) + geom_line() +
  scale_x_continuous(breaks = labels$x, labels = labels$label) +
  ylab("Maximum Take-off Mass (MTOM)") +
  xlab("Projected Axis") + 
  coord_cartesian(ylim = c(5500, 9000)) +
  ggtitle("Effect of Aspect Ratio on MTOM")



## Effect of Lfuse ======================================================================
varL_fusein = seq(10, 15, 0.1)

projectedgraphs <- ProjectedGraphs(xaxis, inp, L_fusein = varL_fusein)
projectedgraphs <- mutate(projectedgraphs, L_fusein = factor(L_fusein, levels = varL_fusein, ordered = TRUE))

labels <- select(projectedgraphs, x, WS, PW) %>% distinct(.keep_all = TRUE) %>% arrange(x) %>%
  mutate(label = paste(round(WS,1), round(PW,1), sep="\n"))
labels <- rbind(data.frame(x = 0, WS = NA, PW = NA, label = "WS:\nPW:"), labels)

ggplot(projectedgraphs, aes(x = x, y = MTOM, group = L_fusein, colour = L_fusein)) + geom_line() +
  scale_x_continuous(breaks = labels$x, labels = labels$label) +
  ylab("Maximum Take-off Mass (MTOM)") +
  xlab("Projected Axis") + 
  coord_cartesian(ylim = c(6800, 9000)) +
  ggtitle("Effect of Length of Fuselage on MTOM")

ggplot(projectedgraphs, aes(x = x, y = Cd0, group = L_fusein, colour = L_fusein)) + geom_line() +
  scale_x_continuous(breaks = labels$x, labels = labels$label) +
  ylab("Maximum Take-off Mass (MTOM)") +
  xlab("Projected Axis") + 
  # coord_cartesian(ylim = c(6800, 9000)) +
  ggtitle("Effect of Length of Fuselage on Cd0")


