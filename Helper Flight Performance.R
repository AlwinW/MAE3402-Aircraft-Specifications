#----------------------------
#--- Flight Performance
#============================

# distance, speed, drag, acceleration, powerout, eta, shaftpower, thrust, lift, drag, Cd0, lift/drag, theta,
# facet 3; start, middle, end; use freex for good view

# Segment | Vinf | Distance | Height | accel | theta_deg | Cl | Cd0 | Cd | L | D | L_D | etaprop | etatotal | PR | PshaftR | TR

selectnames = c("segment", "Vinf", "Distance", "Height", "accel", "theta_deg", "Cl", "Cd0", "L", "D", "L_D", "etaprop", "etatotal", "PR", "PshaftR", "TR")

takeoffdata <- TakeOffLength(inp, all = TRUE)
landingdata <- LandingLength(inp, all = TRUE)
TOall <- filter(takeoffdata$TO, type == "All Engines")
AirDistanceall <- filter(takeoffdata$AirDistanceTO, type == "All Engines")

## Takeoff ======================================================================
varV <- seq(0.1e-5, 1.1*inp$VsTO, length.out = 21)
TOgr <- RepeatRows(inp, varV)
TOgr$varV <- varV

TOgr <- TOgr %>%
  mutate(segment = "Takeoff") %>%
  rowwise() %>%
  do(data.frame(
    .,
    Distance = integrate(function(V) 
      V/GroundAcceleration(TOall, V, distancecalc = TRUE), 
      lower = 0, .$varV)[[1]]
  )) %>%
  rowwise() %>%
  do(data.frame(
    .,
    GroundAcceleration(TOall, .$varV, all = TRUE)
  )) %>%
  mutate(Vinf = velval,
         Height = 0,
         theta_deg = 0,
         qinf = 1/2*rho*Vinf^2,
         L = qinf * S * Cl,
         D = qinf * S * Cd,
         L_D = L/D,
         etatotal = etamech * etaprop,
         PR = PA,
         PshaftR = PA/etaprop,
         TR = TA)
TOgr <- data.frame(TOgr[selectnames])
# You can do freeroll properly later if you want
TOgr <- rbind(
  TOgr,
  mutate(TOgr[nrow(TOgr),], 
         Distance = TOgr[nrow(TOgr),"Distance"] + 1.1 * 3 * inp$VsTO)
  )


varh <- seq(0, AirDistanceall$hTR, length.out = 11)
TOtr <- RepeatRows(inp, varh)
TOtr$h <- varh
TOtr$gamma <- AirDistanceall$gamma

TOtr <- StandardAtomsphere(TOtr) %>%
  mutate(
    segment = "Takeoff",
    Height = h,
    Vinf = 1.15 * VsTO,
    qinf = 1/2 * rho * Vinf^2,
    Cd0 = Cd0clean + Cd0flaps + Cdiflaps,
    theta_deg = gamma *180/pi,
    L = W/cos(gamma),
    Cl = L / (qinf * S),
    Cd = Cd0 + K * Cl^2,
    D = qinf * S * Cd,
    L_D = L/D,
    Distance = h/tan(gamma),
    accel = 0,
    etaprop = etaprop(Vinf),
    etatotal = etamech * etaprop,
    PR = D + W*sin(gamma),
    PshaftR = PR/etaprop,
    TR = PR/Vinf
  )
TOtr <- data.frame(TOtr[selectnames])


## Segments ======================================================================
varh <- seq(AirDistanceall$hTR, inp$AltFlaps,length.out = 11)
Seg2 <- RepeatRows(inp, varh)
Seg2$h <- varh
Seg2$gamma <- AirDistanceall$gamma
Segment2 <- 

