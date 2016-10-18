#----------------------------
#--- Energy Usage
#============================

ClimbEnergy <- function(inp, h, energycalc = FALSE, sinonly = FALSE) {
  inph <- RepeatRows(inp, length(h))
  inph$h <- h
  inph <- StandardAtomsphere(inph) %>%
    mutate(
      Vinf = ifelse(Vconst != 0, Vconst, Mach*a),
      etaprop = etaprop(Vinf),
      PA = PA(Pshafteng, Ne, Vinf)) %>%
    rowwise() %>%
    do(data.frame(
      .,
      ClimbRatesFunction(.$PA, .$Cd0, .$rho, .$Vinf, .$S, .$K, .$W)
    )) %>%
    ungroup() %>%
    mutate(
      qinf = 1/2 * rho * Vinf^2,
      Cl = nload *W/(qinf * S),
      Cd = Cd0 + K * Cl^2,
      D = Cd * qinf * S,
      TR = D + W * SinTheta,
      PR = TR * Vinf,
      P_Vv = (PR/etaprop)/ClimbRate
    )
  data.frame(inph)
  if (energycalc == TRUE)
    return(inph$P_Vv)
  else if (sinonly == TRUE)
    return(inph$SinTheta)
}

AccelerateEnergy <- function(inp, V, energycalc = FALSE) {
  inpv <- RepeatRows(inp, length(V))
  inpv$Vinf<- V
  inpv <- inpv %>%
    mutate(
      qinf = 1/2 * rho * Vinf^2,
      Cl = W/(qinf * S),
      Cd = Cd0 + K*Cl^2,
      D = Cd*qinf * S,
      accel = D/m,
      Dv_accel = D * Vinf/accel
    )
  if (energycalc == TRUE)
    return(inpv$Dv_accel)
}

DescentRatesFunction <- function(Cd0, K, W, m, S, h, V, theta, 
                                 accelonly = FALSE, all = FALSE, Vvonly = FALSE) {
  rho = StandardAtomsphere(data.frame(h=h))$rho
  qinf = 1/2 * rho * V^2
  Cl = W*cos(theta)/(qinf*S)
  Cd = Cd0 + K * Cl^2
  D = qinf*Cd*S
  accel = (W*sin(theta) - D)/m
  Vv = -V*sin(theta)
  Vh = V*cos(theta)
  accelv = -accel*sin(theta)
  accelh = accel*cos(theta)
  if (accelonly == TRUE)
    return(accel)
  if (Vvonly == TRUE)
    return(Vv)
  if (all == TRUE)
    return(data.frame(
      rho = rho,
      Cl = Cl,
      Cd = Cd,
      D = D,
      accel = accel,
      Vv = Vv,
      Vh = Vh,
      accelv = accelv,
      accelh = accelh
    ))
}

DescentZeroAccel <- function(inp, h, V, timecalc = FALSE, all = TRUE) {
  Desc <- RepeatRows(inp, length(h))
  Desc$h <- h
  Desc$V <- V
  # Find theta for zero accel
  Desc <- Desc %>%
    rowwise() %>%
    do(data.frame(
      .,
      theta = ModifiedSecant(function(theta) 
        DescentRatesFunction(.$Cd0, .$K, .$W, .$m, .$S, .$h, .$V, theta, accelonly = TRUE),
        xr = 0.063, del = 1e-10, toler = 1e-8, positive = TRUE)
    )) %>%
    rowwise() %>%
    do(data.frame(
      .,
      DescentRatesFunction(.$Cd0, .$K, .$W, .$m, .$S, .$h, .$V, .$theta, all = TRUE)
    )) %>%
    mutate(accel = 0)
  Desc <- data.frame(Desc)
  if (timecalc == TRUE)
    return(1/Desc$Vv)
  else if (all == TRUE)
  return(Desc)
}

## Energy Summary ======================================================================
EnergyUsage <- function(TO, AirDistance, inp) {
## Take Off ======================================================================
  #--- Clean up input
  TOall <- filter(TO, type == "All Engines")
  AirDistanceall <- filter(AirDistance, type == "All Engines")
  #--- Take-off Ground Roll
  TOgr <- integrate(function(V)
      GroundAcceleration(TOall, V, energycalc = TRUE),
    lower = 0, upper = 1.1*inp$VsTO
  )[[1]] +
    PA(inp$Pshafteng, 2, 1.1 * inp$VsTO) * 3
  #--- Take-off Transition
  TOtr <- integrate(Vectorize(function(h)
    PA(inp$Pshafteng, 2, 1.15*inp$VsTO) / 
      (sin(AirDistanceall$gamma) * 1.15*inp$VsTO)),
    lower = 0, upper = AirDistance$hTR
  )[[1]]
## Climb Segments ======================================================================
  #--- 1st Semgnet Climb
  if (AirDistanceall$SC > 0) {
    Seg1 <- integrate(Vectorize(function(h)
      PA(inp$Pshafteng, 2, 1.15*inp$VsTO) / 
        (sin(atan(AirDistanceall$hTR/(AirDistanceall$Sair - AirDistanceall$ST))) * 1.15*inp$VsTO)),
      lower = AirDistance$hTR, upper = inp$Hobs
    )[[1]]
  } else {
    Seg1 = 0
  }
  #--- 2nd Segment Climb
  Seg2inp <- inp %>%
    mutate(
      Cd0 = Cd0clean + Cd0flaps + Cdiflaps,
      Ne = 2,
      Vconst = 1.2*VsTO
    )
  Seg2 <- integrate(function(h)
    ClimbEnergy(Seg2inp, h, energycalc = TRUE),
    lower =  max(AirDistance$hTR, inp$Hobs), upper = inp$AltFlaps
    )[[1]]
  #--- 3rd Segment Acceleration
  Seg3inp <- inp %>%
    mutate(
      h = AltFlaps,
      Cd0 = Cd0clean + Cd0flaps + Cdiflaps,
      Ne = 2
    )  %>%
    StandardAtomsphere(.)
  Seg3 <- integrate(function(V)
    AccelerateEnergy(Seg3inp, V, energycalc = TRUE),
    lower = 1.2*inp$VsTO, upper = inp$Mach*Seg3inp$a
  )[[1]]
  # 4th Segment Climb
  Seg4inp <- inp %>%
    mutate(
      Cd0 = Cd0clean,
      Ne = 2,
      Vconst = 0
    )
  Seg4 <- integrate(function(h)
    ClimbEnergy(Seg4inp, h, energycalc = TRUE),
    lower =  inp$AltFlaps, upper = inp$AltCruise
  )[[1]]
## Descent Segments ======================================================================
  #--- Critical h Value
  Desc <- mutate(inp, Cd0 = Cd0clean + 2*Cd0propfea)
  hcrit = ModifiedSecant(function(h) 
    DescentRatesFunction(Desc$Cd0, Desc$K, Desc$W, Desc$m, Desc$S, h, Desc$Vcruise, Desc$DescAngle, accelonly = TRUE),
    xr = 1000, del = 0.001, toler = 0.01, positive = TRUE)
  hcrit = max(hcrit, inp$AltFlaps)
  #--- Zero Acceleration
  timecrit = integrate(function(h)
    DescentZeroAccel(Desc, h, inp$Vcruise, timecalc = TRUE),
    lower = inp$AltCruise, upper = hcrit
    )[[1]]
  timeflap = integrate(function(h)
    DescentRatesFunction(Desc$Cd0, Desc$K, Desc$W, Desc$m, Desc$S, h, Desc$Vcruise, Desc$DescAngle, accelonly = TRUE),
    lower = hcrit, upper = inp$AltFlaps
  )[[1]]
  
  time = integrate(function(h)
    DescentZeroAccel(Desc, h, inp$Vcruise, timecalc = TRUE),
    lower = inp$AltCruise, upper = inp$AltFlaps
  )[[1]]
  
  
  
  
  energysummary <- data.frame(
    TOgr = TOgr,
    TOtr = TOtr,
    Seg1 = Seg1,
    Seg2 = Seg2,
    Seg3 = Seg3,
    Seg4 = Seg4
  )
  return(energysummary)
}