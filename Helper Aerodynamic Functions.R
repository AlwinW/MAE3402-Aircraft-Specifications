#----------------------------
#--- Aerodynamic Calculations
#============================

## Calculate Cd given Cl ======================================================================
Cd <- function(Cd0, K, Cl)
  Cd0 + K * Cl^2

## Min Velocity for Cl ======================================================================
Vmin <- function(rho, WS, Clmax) 
  sqrt(2/rho * WS * 1/Clmax)

## Power Derating Function ======================================================================
PA <- function(Pshafteng, Ne, Vinf)
  Pshafteng * Ne * etaprop(Vinf)

## K Effective Due to Ground Effect ======================================================================
Keff <- function(K, h, b)
  (33 * (h/b)^1.5) / (1 + 33 * (h/b)^1.5)  * K

## Estimated Propeller Efficiency ======================================================================
etaprop <- function(V)
  -0.000117737 * V^2 + 0.02054599 * V

## Friction on the Ground ======================================================================
groundmu <- data.frame(names = c("Dry Concrete", "Wet Concrete", "Icy Concrete"),
                       brakesoff = c(0.04, 0.05, 0.02), #NB: 0.03-0.05 for dry
                       brakeson_min = c(0.3, 0.15, 0.06),
                       brakeson_max = c(0.5, 0.3, 0.10)) %>%
  mutate(brakeson = (brakeson_min + brakeson_max)/2)
rownames(groundmu) <- groundmu$names

## Climb Rates for Given Power,  etc ======================================================================
ClimbRatesFunction <- function(P, Cd0, rho, V, S, K, W) {
  # Coefficients to simplify the calculation
  a = P - 1/2 * Cd0 * rho * V^3 * S
  b = (2 * K * W^2) / (rho * V * S)
  c = W * V
  # Solution to the quadratic
  sintheta = c((c - sqrt(-4*a*b + 4*b^2 + c^2)) / (2*b),
               (c + sqrt(-4*a*b + 4*b^2 + c^2)) / (2*b))
  # Determine which angle to accept
  sintheta = sintheta[sintheta < 0.5 & sintheta > -0.5]
  # Return a meaningful result
  if (length(sintheta) != 1) {
    return(data.frame(Theta = NA, SinTheta = NA, PerGrad = NA, ClimbRate = NA))
  }
  else {
    theta = asin(sintheta)
    return(data.frame(
      Thetadeg = theta * 180 / pi,  # NOTE: Thetadeg is in degrees!! theta is in radians
      SinTheta = sintheta, 
      PerGrad = tan(theta)*100, 
      ClimbRate = sintheta * V,
      nload = 1/cos(theta) * sign(theta)
      ))
  }
}

## Aerodyanmic Parameters at Various Segments ======================================================================
SegmentParameters <- function(inp) {
  # Use recode {dplyr} if you have time...
}

ClimbEnergy <- function(inp, h, energycalc = FALSE, sinonly = FALSE, distancecalc = FALSE) {
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
      P_Vv = (PR/etaprop)/ClimbRate,
      accel = TR/m,
      Vh = Vinf*cos(Thetadeg * pi/180),
      distcalc = Vh/ClimbRate
    )
  data.frame(inph)
  if (energycalc == TRUE)
    return(inph$P_Vv)
  else if (sinonly == TRUE)
    return(inph$SinTheta)
  else if (distancecalc == TRUE)
    return(inph$distcalc)
}

AccelerateEnergy <- function(inp, V, energycalc = FALSE, distancecalc = FALSE) {
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
  if (distancecalc == TRUE)
    return(inpv$accel)
}

DescentRatesFunction <- function(Cd0, K, W, m, S, h, V, theta, 
                                 distancecalc = FALSE,
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
  distcalc = Vh/Vv
  if (distancecalc == TRUE)
    return(distcalc)
  else if (accelonly == TRUE)
    return(accel)
  else if (Vvonly == TRUE)
    return(Vv)
  else if (all == TRUE)
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

DescentZeroAccel <- function(inp, h, V, timecalc = FALSE, distancecalc = FALSE, all = FALSE) {
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
    mutate(accel = 0,
           Vh = V*cos(theta * pi/180),
           distcalc = Vh/Vv)
  Desc <- data.frame(Desc)
  if (timecalc == TRUE)
    return(1/Desc$Vv)
  else if (distancecalc == TRUE)
    return(Desc$distcalc)
  else if (all == TRUE)
    return(Desc)
}
