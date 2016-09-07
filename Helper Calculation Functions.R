#----------------------------
#--- Calculations Functions
#============================

## Calculate Cd given Cl ======================================================================
Cd <- function(Cd0, K, Cl)
  Cd0 + K * Cl^2 # Parabolic Relation

## Min Velocity for Cl ======================================================================
Vmin <- function(rho, WS, Clmax) 
  sqrt(2/rho * WS * 1/Clmax)

## Power Derating Function ======================================================================
PA <- function(P0eng, sigma)
  P0eng * sigma^inputvals$alt_s

## K Effective Due to Ground Effect ======================================================================
Keff <- function(K, h, b)
  (33 * (h/b)^1.5) / (1 + 33 * (h/b)^1.5)  * K

## Friction on the Ground ======================================================================
groundmu <- data.frame(names = c("Dry Concrete", "Wet Concrete", "Icy Concrete"),
                       brakesoff = c(0.04, 0.05, 0.02), #NB: 0.03-0.05 for dry
                       brakeson_min = c(0.3, 0.15, 0.06),
                       brakeson_max = c(0.5, 0.3, 0.10)) %>%
  mutate(brakeson = (brakeson_min + brakeson_max)/2)

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
      Theta = theta * 180 / pi, 
      SinTheta = sintheta, 
      PerGrad = tan(theta)*100, 
      ClimbRate = sintheta * V))
  }
}

## Ground Roll Distances ======================================================================
# Accelerate to V1 then brake to stop
AccelerateStop <- function(coef, AirDistance, V1, V2) {
  integrate(function(x) x^2 / (coef$A[1] + coef$B[1] * x + coef$C[1] * x^3), 0, V1)[[1]] +
    V1 +
    integrate(function(x) x^2 / (coef$A[3] + coef$B[3] * x + coef$C[3] * x^3), V1, 0)[[1]]
  }
# Accelerate to V1 then continue after engine failure to V2 + 3 seconds reaction + air distance 
AccelerateContinue <- function(coef, AirDistance, V1, V2) {
  integrate(function(x) x^2 / (coef$A[1] + coef$B[1] * x + coef$C[1] * x^3), 0, V1)[[1]] +
    integrate(function(x) x^2 / (coef$A[2] + coef$B[2] * x + coef$C[2] * x^3), V1, V2)[[1]] +
    3 * V2 +
    AirDistance$Sair[2]
}
# Accelerate to V2 + 3 seconds reaction + air distance
AccelerateLiftOff <- function(coef, AirDistance, V1, V2) {
  (integrate(function(x) x^2 / (coef$A[1] + coef$B[1] * x + coef$C[1] * x^3), 0, V2)[[1]] +
     3 * V2 +
     AirDistance$Sair[1]) * 1.15
}
# Deccelerate from Vtd to 0 + 3 seconds reaction + air distance
DeccelerateStop <- function(coef, AirDistance, Vtd) {
  (integrate(function(x) x^2 / (coef$A[1] + coef$B[1] * x + coef$C[1] * x^3), Vtd, 0)[[1]] +
     3 * Vtd +
     AirDistance$Sair) * 1.67
}

## Main: Combine Input and Specs ======================================================================
InpSpecs <- function(inputvals, specifications) {
  inp  <- t(specifications["Value"])
  colnames(inp) <- t(specifications["Variable"])
  inp <- cbind(inputvals, inp)
  return(inp)
}

## Main: AeroParams ======================================================================
AeroParams <- function(inp, iteration = FALSE) {
  out <-  inp[rep(row.names(inp), each = 5), 1:length(inp)]
  out$type <- c("Sea Level", "Cruise", "Ceiling", "Takeoff", "Landing")
  out$h <- c(0, inp$AltCruise, inp$AltCeil, 0, 0)
  out$Clmax <- inp$Clclean + c(0, 0, 0, inp$Clflaps, inp$Clhls)
  # Transform the dataframe to find the various aerodynamic properties
  out <- StandardAtomsphere(out) %>%
    mutate(Vinf = Mach * a,
           Vstall = Vmin(rho, WS, Clmax))
  out$Vinf <- c(out$Vinf[1], out$Vinf[2], out$Vinf[3], out$Vstall[4] * 1.2, out$Vstall[5] * 1.3)
  out <- out %>%
    mutate(Vsafe = 1.2 * Vstall,
           qinf = 1/2 * rho * Vinf^2,
           Cl = W / (qinf * S),
           Cd = Cd0 + K * Cl^2,
           ClCd = Cl/Cd,
           Clstar = sqrt(Cd0 / K),
           Cdstar = 2 * Cd0,
           ClCdstar = 1 / sqrt(4 * Cd0 * K),
           Vstar = Vmin(rho, WS, Clstar),
           Cl32 = sqrt(3) * Clstar,
           Cd32 = 2 * Cdstar,
           ClCd32 = sqrt(3/4) * ClCdstar,
           V32 = (1/3)^(1/4) * Vstar
    )
  #--- Create a meaningful output to be returned  to a user
  AeroParamsTable <- select(
    out,
    type, h, rho, Vinf, Vstall, Vsafe, Vstar, V32, 
    Cl, Clstar, Cl32, Clmax, Cd, Cdstar, Cd32, 
    ClCd, ClCdstar, ClCd32)
  #--- Return the result as a list
  return(list(out = out, AeroParamsTable = AeroParamsTable))
}

## Main: Takeoff ======================================================================
TakeOff <- function(inp, iteration = FALSE) {
  #--- Initialise a data frame to apply functions to
  coef <- inp[rep(row.names(inp), each = 3), 1:length(inp)]
  coef$type <- c("All Engines", "One Engine Down", "Rejected Take-Off")
  coef$Ne <- c(2, 1, 0)
  coef$mu <- c(
    as.double(filter(groundmu,names == "Dry Concrete") %>% select(brakesoff)),
    as.double(filter(groundmu,names == "Dry Concrete") %>% select(brakesoff)),
    as.double(filter(groundmu,names == "Dry Concrete") %>% select(brakeson))
  )
  #--- Determine the important coefficients and values
  coef <- mutate(coef, h = 0) %>%
    StandardAtomsphere(.) %>%
    mutate(
      Cl = ClG,
      Keff = Keff(K, hground, b),
      Cd  = Cd0G + Keff * ClG ^ 2,
      Vstall = Vmin(rho, WS, Clclean + Clflaps),
      Vlof = 1.1 * Vstall,
      A = g * (PA(P0eng, sigma) * Ne / W),
      B = g * (- mu),
      C = g * (rho/2 * 1/WS * (mu * Cl - Cd))
    )
  V2 <- coef$Vlof[1]
  #--- Determine the AirDistance in each case
  AirDistance <- inp[rep(row.names(inp), each=2), 1:length(inp)]
  rownames(AirDistance) <- NULL
  AirDistance$type <- c("All Engines", "One Engine Down")
  AirDistance$Ne <- c(2, 1)
  AirDistance <- AirDistance %>%
    mutate(h = 0) %>%
    StandardAtomsphere(.) %>%
    mutate(ClTR = Clclean + Clflaps, 
           Vstall = Vmin(rho, WS, ClTR), VTR = Vstall * 1.15,
           PA = PA(P0eng * Ne, sigma), TA = PA / VTR,
           qinf = 1/2 * rho * VTR^2, Cd = Cd0 + K * ClTR^2,
           D = qinf * S * Cd, L = qinf * S * ClTR) %>%
    rowwise() %>%
    mutate(
      R = (VTR) ^ 2 / (0.2 * g),
      gamma = ClimbRatesFunction(PA, Cd0, rho, VTR, S, K, W)[[1]] * pi / 180,
      hTR = R * (1 - cos(gamma)),
      ST = R * (sin(gamma)),
      SC = (Hobs - hTR) / tan(gamma),
      Sair = ifelse(SC >0, ST + SC, sqrt(R^2 - (R-hTR)^2))
    )%>%
    ungroup()
  AirDistance <- data.frame(select(AirDistance, type, R, gamma, hTR, ST, SC, Sair))
  #--- Determine the BFL
  ## NOTE I SHOULD ALLOW FOR 1s FOR THE PILOT TO REALISE THE ENGINE FAILURE (RAYMER)
  V1 = ModifiedSecant(function(x) AccelerateStop(coef, AirDistance, x, V2) - AccelerateContinue(coef, AirDistance, x, V2), 
                      V2, 0.01, 1e-4, positive = TRUE)
  #--- Output the data
  return(
    data.frame(
      NTO = AccelerateLiftOff(coef, AirDistance, V1, V2),
      BFL = AccelerateStop(coef, AirDistance, V1, V2),
      V1 = V1,
      V2 = V2) %>%
      mutate(TakeOffDistance = max(NTO, BFL))
  )
}

## Main:Climb ======================================================================
#--- Determine the various climb rates required
ClimbRates <- function(inp, iteration = FALSE) {
  # Create a data frame of the three scenarios
  out <-  inp[rep(row.names(inp), each = 3), 1:length(inp)]
  out$type <- c("2nd Seg OEI Climb", "Cruise", "Ceiling")
  out$Ne <- c(1, 2, 2)
  out$h <- c(inp$Hobs, inp$AltCruise, inp$AltCeil)
  out$Clmax <- inp$Clclean + c(inp$Clflaps, 0, 0)
  # Determine the climb rates for each scenario
  out <- StandardAtomsphere(out) %>%
    mutate(Vinf = Mach * a,
           Vstall = Vmin(rho, WS, Clmax),
           Vsafe = 1.2 * Vstall)
  out$Vinf <- c(out$Vsafe[1], out$Vinf[1], out$Vinf[1])
  out <- mutate(
      out,
      qinf = 1/2 * rho * Vinf^2,
      Cl = W / (qinf * S),
      Cd = Cd0 + K * Cl^2,
      PA = PA(P0eng, sigma) * Ne) %>%
    rowwise() %>%
    do(data.frame(., ClimbRatesFunction(.$PA, .$Cd0, .$rho, .$Vinf, .$S, .$K, .$W))) %>%
    ungroup()
  return(data.frame(out))
}

## Main: Landing ======================================================================
#--- Determine the distance for landing
Landing <- function(inp, iteration = FALSE) {
  #--- Determine the distance required for landing
  AirDistance <- inp
  AirDistance$type <- "All Engines"
  AirDistance$Ne <- 2
  AirDistance <- AirDistance %>%
    mutate(h = 50*0.3048) %>%
    StandardAtomsphere(.) %>%
    mutate(
      Clmax = Clclean + Clhls,
      Vstall = Vmin(rho, WS, Clmax),
      Vapp = 1.3* Vstall,
      qinf = 1/2 * rho * Vapp^2,
      gamma = max(-3, ClimbRatesFunction(0.01*P0, Cd0G, rho, Vapp, S, K, W)[[1]]),
      L = W * cos(gamma * pi / 180),
      Cl = L / (qinf * S),
      Cd = Cd0G + K * Cl^2,
      D = qinf * S * Cd,
      TR = D - W * sin(gamma * pi / 180),
      PR = TR * Vapp,
      R = Vapp^2 / (0.2 * g),
      SF = R * sin(gamma * pi / 180),
      hF = R * (1 - cos(gamma * pi / 180)),
      SA = (50*0.3048 - hF) / tan(gamma*pi/180),
      Sair = ifelse(SA >0, SA + SF, sqrt(R^2 - (R-50*0.3048)^2)))

  #--- Initialise a data frame to apply functions to
  Vapp = AirDistance$Vapp # N.B @ h = 50ft
  coef <- inp
  coef$type <- "Engines off"
  coef$Ne <- 0
  coef$mu <- as.double(filter(groundmu, names == "Dry Concrete") %>% select(brakeson))
  coef <- mutate(coef, h = 0) %>%
    StandardAtomsphere(.) %>%
    mutate(Vstall = Vmin(rho, WS, Clclean + Clhls),
           Vtd = 1.15 * Vstall,
           Vapp = 1.3 * Vstall)
  Vtd = coef$Vtd
  #--- Determine the important coefficients and values
  coef <- coef %>%
    mutate(
      Cl = ClG,
      Keff = Keff(K, hground, b),
      Cd  = Cd0G + Keff * ClG ^ 2,
      A = g * (PA(P0eng, sigma) * Ne / W),
      B = g * (- mu),
      C = g * (rho/2 * 1/WS * (mu * Cl - Cd))
    )
  #--- Output the data
  return(
    data.frame(
      Vtd = Vtd,
      Vapp = Vapp,
      LandingDistance = DeccelerateStop(coef, AirDistance, Vtd))
  )
  
}

## Main: Power Estimate ======================================================================
#--- Estimate the power usage with Range only
PowerEst <- function(inp, iteration = FALSE) {
  out <- mutate(inp, h = AltCruise) %>%
    StandardAtomsphere(.) %>%
    mutate(Vinf = Mach * a,
           qinf = 1/2 * rho * Vinf^2,
           Cl = W/(qinf * S),
           Cd = Cd0 + K*Cl^2,
           D = qinf * Cd * S,
           Wb = D/Etatotal * 1.05,
           Vb = Wb / Dens,
           WbW0 = Wb/m)
  Wbest <- select(out, Wb, Vb, WbW0)
  return(list(out = out, Wbest = Wbest))
}

## Main: Power Values ======================================================================
#--- Determine the power usage in each section
PowerUsage <- function(inp, iteration = FALSE, resolution = 10) {
  #--- Summary of takeoff
  Takeoffout <- TakeOff(inp)
  #--- Initialise the takeoff
  PTOgr <- mutate(inp, h = 0) %>%
    StandardAtomsphere(.) %>%
    mutate(
      p = 1/2 * rho * S,
      Cl = ClG,
      Keff = Keff(K, hground, b),
      Cd  = Cd0G + Keff * ClG ^ 2)
  PTOgr$Ne <- 2
  PTOgr$mu <- as.double(filter(groundmu,names == "Dry Concrete") %>% select(brakesoff))

  
  LevelRollInt <- function(inp, V){
    Int <- inp[rep(row.names(inp), each = length(V)), 1:length(inp)]
    Int$V <- V
    Int = mutate(Int,
      TR = p*V^2*Cd + mu*(W - p*V^2*Cl),
      PR = TR * V,
      accel = (PA(P0eng, sigma) * Ne/V - p*V^2*Cd - mu*(W - p*V^2*Cl)) / m,
      Int =  PR/ accel)
    return(Int)
  }
  
  integrate(function(x) LevelRollInt(PTOgr, x)$Int, 0, Takeoffout$V2)[[1]]
  
  asdf <- LevelRollInt(PTOgr, seq(1e-1,Takeoffout$V2, length.out = resolution * 5))
  
  asdf <- asdf %>%
    rowwise() %>%
    mutate(
      Distance = integrate(function(x) LevelRollInt(PTOgr, x)$V / LevelRollInt(PTOgr, x)$accel, 0, V)[[1]],
      Time = integrate(function(x) 1/LevelRollInt(PTOgr, x)$accel, 10e-8, V)[[1]],
      Energy = integrate(function(x) LevelRollInt(PTOgr, x)$PR / LevelRollInt(PTOgr, x)$accel, 0, V)[[1]]
      )
  asdf <- data.frame(asdf)
  
}


