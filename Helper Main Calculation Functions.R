#----------------------------
#--- Main Calculation Functions
#============================

## Groundroll ======================================================================
GroundAcceleration <- function(inp, velval, intval = TRUE) {
  inp <- RepeatRows(inp, velval)
  inp$velval = velval
  inp <- inp %>%
    mutate(
      qinf = 1/2 * rho * velval^2,
      D = qinf * Cd * S,
      L = qinf * Cl * S,
      Fr = mu * (W - L),
      etaprop = etaprop(velval),
      PA = PA(Pshafteng, Ne, velval),
      TA = PA / velval,
      Fnet = (TA - D - Fr),
      accel = Fnet/m
    )
  if (intval == TRUE)
    return(inp$accel)
  else
    return(inp)
}

## Takeoff Helper Functions ======================================================================
# Accelerate to V1 then brake to stop
AccelerateStop <- function(TO, V1) {
  integrate(function(V) V/GroundAcceleration(filter(TO, type == "All Engines"), V), 0, V1)[[1]] + 
    V1 * (1) +
    integrate(function(V) V/GroundAcceleration(filter(TO, type == "Rejected Take-Off"), V), V1, 0)[[1]]
}
# Accelerate to V1 then continue after engine failure to V2 + 3 seconds reaction + air distance 
AccelerateStop <- function(TO, AirDistance, V1, V2) {
  integrate(function(V) V/GroundAcceleration(filter(TO, type == "All Engines"), V), 0, V1)[[1]] + 
    integrate(function(V) V/GroundAcceleration(filter(TO, type == "One Engine Down"), V), V1, V2)[[1]] + 
    V2 * (3) +
    AirDistance$Sair[2]
}

## Takeoff ======================================================================
# TakeOff <- function(inp) {
  #--- Set up the initial parameters to solve for
  TO <- RepeatRows(inp, 3)
  TO$segment <- "Takeoff"
  TO$type <- c("All Engines", "One Engine Down", "Rejected Take-Off")
  TO$Ne <- c(2, 1, 0)
  TO$mu <- c(
    as.double(groundmu["Dry Concrete", "brakesoff"]),
    as.double(groundmu["Dry Concrete", "brakesoff"]),
    as.double(groundmu["Dry Concrete", "brakeson"])
  )
  #---- Determine aerodynamic parameters
  TO <- mutate(TO, h = 0) %>%
    StandardAtomsphere(.) %>%
    mutate(
      Cd0 = Cd0clean + Cd0lg + Cd0flaps + Cdiflaps,
      Cd0 = Cd0 + as.numeric(Ne == 1) * Cd0propfea + as.numeric(Ne == 0) * 2 * Cd0propunfea,
      Cl = Cl0 + Clflaps
    )
  #--- Determine the AirDistance in each case
  AirDistance <- TO %>%
    filter(type %in% c("All Engines", "One Engine Down")) %>%
    mutate(ClTR = Clclean + Clflaps, 
           VTR = VsTO * 1.15,
           etaprop = etaprop(VTR),
           PA = PA(Pshafteng, Ne, VTR), 
           TA = PA / VTR,
           qinf = 1/2 * rho * VTR^2, 
           Cd = Cd0 + K * ClTR^2,
           D = qinf * S * Cd, 
           L = qinf * S * ClTR) %>%
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
  #--- Calculate the other TO parameters and find V1
  TO <- TO %>%
    mutate(
      K = Keff(K, hground, b),
      Cd = Cd0 + K * Cl^2
    )
# }