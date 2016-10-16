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
AccelerateContinue <- function(TO, AirDistance, V1, V2) {
  integrate(function(V) V/GroundAcceleration(filter(TO, type == "All Engines"), V), 0, V1)[[1]] + 
    integrate(function(V) V/GroundAcceleration(filter(TO, type == "One Engine Down"), V), V1, V2)[[1]] + 
    V2 * (3) +
    AirDistance$Sair[2]
}
# Accelerate to V2 + 3 seconds reaction + air distance * Safety Factor 
AccelerateLiftOff <- function(TO, AirDistance, V2) {
  (integrate(function(V) V/GroundAcceleration(filter(TO, type == "All Engines"), V), 0, V2)[[1]] +
     3 * V2 +
     AirDistance$Sair[1]) * 1.15
}

## Takeoff ======================================================================
TakeOffLength <- function(inp, V2 = 1.1 * inp$VsTO) {
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
  #--- Test if the graphs will cross over before V2
  BFL <- try(ModifiedSecant(function(V1) AccelerateContinue(TO, AirDistance, V1, V2) - AccelerateStop(TO, V1), 
                             V2, 0.01, 1e-4, positive = TRUE), silent = TRUE)
  if (is.numeric(BFL)){
    V1 = ModifiedSecant(function(V1) AccelerateContinue(TO, AirDistance, V1, V2) - AccelerateStop(TO, V1), 
                   V2, 0.01, 1e-4, positive = TRUE)
    BFL = AccelerateContinue(TO, AirDistance, V1, V2)
  } else {
    V1 = NA
    BFL = 0
  }
  #--- Output the results
  TOoutput <- data.frame(
    V1 = V1,
    V2 = V2,
    BFL = BFL,
    NTO = AccelerateLiftOff(TO, AirDistance, V2)
  ) %>%
    mutate(
      NTOgr = NTO - AirDistance$Sair[1],
      TakeOffLength = max(BFL*as.numeric(V1 <= V2), NTO)
    )
  return(TOoutput)
}
#--- Plot of the takeoff curves
TakeOffLengthPlot<- function(TO, AirDistance, TOoutput, inp) {
  TOplotdata <- data.frame(
    V1 = seq(0.00001, max(TOoutput$V1, TOoutput$V2)*1.02, len = 20),
    V2 = TOoutput$V2) %>%
    rowwise() %>%
    do(data.frame(
      .,
      AccelerateStop = AccelerateStop(TO, .$V1),
      AccelerateContinue = AccelerateContinue(TO, AirDistance, .$V1, .$V2)
      ))
  TOplotout <- ggplot(data = TOplotdata, aes(x = V1)) +
    geom_hline(yintercept = inp$Srun, aes(colour = "Max Paved Runway")) + 
    geom_line(aes(y = AccelerateStop, colour = "Accelerate-Stop")) + 
    geom_line(aes(y = AccelerateContinue, colour = "Accelerate-Continue"))
}

## Climb ======================================================================











