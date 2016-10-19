#----------------------------
#--- Main Calculation Functions
#============================

## Groundroll ======================================================================
GroundAcceleration <- function(TO, velval, distancecalc = FALSE, energycalc = FALSE, all = FALSE) {
  TOgr <- RepeatRows(TO, length(velval))
  TOgr$velval = velval
  TOgr <- TOgr %>%
    mutate(
      qinf = 1/2 * rho * velval^2,
      D = qinf * Cd * S,
      L = qinf * Cl * S,
      Fr = mu * (W - L),
      etaprop = etaprop(velval),
      PA = PA(Pshafteng, Ne, velval),
      TA = PA / velval,
      Fnet = (TA - D - Fr),
      accel = Fnet/m,
      Pshaft_accel = (PA/etaprop)/accel
    )
  if (distancecalc == TRUE)
    return(TOgr$accel)
  else if (energycalc == TRUE)
    return(TOgr$Pshaft_accel)
  else if (all == TRUE)
    return(TOgr)
}

## Takeoff Helper Functions ======================================================================
# Accelerate to V1 then brake to stop
AccelerateStop <- function(TO, V1) {
  integrate(function(V) V/GroundAcceleration(filter(TO, type == "All Engines"), V, distancecalc = TRUE), 0, V1)[[1]] + 
    V1 * (1) +
    integrate(function(V) V/GroundAcceleration(filter(TO, type == "Rejected Take-Off"), V, distancecalc = TRUE), V1, 0)[[1]]
}
# Accelerate to V1 then continue after engine failure to V2 + 3 seconds reaction + air distance 
AccelerateContinue <- function(TO, AirDistance, V1, V2) {
  integrate(function(V) V/GroundAcceleration(filter(TO, type == "All Engines"), V, distancecalc = TRUE), 0, V1)[[1]] + 
    integrate(function(V) V/GroundAcceleration(filter(TO, type == "One Engine Down"), V, distancecalc = TRUE), V1, V2)[[1]] + 
    V2 * (3) +
    AirDistance$Sair[2]
}
# Accelerate to V2 + 3 seconds reaction + air distance * Safety Factor 
AccelerateLiftOff <- function(TO, AirDistance, V2) {
  (integrate(function(V) V/GroundAcceleration(filter(TO, type == "All Engines"), V, distancecalc = TRUE), 0, V2)[[1]] +
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
  AirDistanceTO <- TO %>%
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
  AirDistanceTO <- data.frame(select(AirDistanceTO, type, R, gamma, hTR, ST, SC, Sair))
  #--- Calculate the other TO parameters and find V1
  TO <- TO %>%
    mutate(
      K = Keff(K, hground, b),
      Cd = Cd0 + K * Cl^2
  )
  #--- Test if the graphs will cross over before V2
  BFL <- try(ModifiedSecant(function(V1) AccelerateContinue(TO, AirDistanceTO, V1, V2) - AccelerateStop(TO, V1), 
                             V2, 0.01, 1e-4, positive = TRUE), silent = TRUE)
  if (is.numeric(BFL)){
    V1 = ModifiedSecant(function(V1) AccelerateContinue(TO, AirDistanceTO, V1, V2) - AccelerateStop(TO, V1), 
                   V2, 0.01, 1e-4, positive = TRUE)
    BFL = AccelerateContinue(TO, AirDistanceTO, V1, V2)
  } else {
    V1 = NA
    BFL = 0
  }
  #--- Output the results
  TOoutput <- data.frame(
    V1 = V1,
    V2 = V2,
    BFL = BFL,
    NTO = AccelerateLiftOff(TO, AirDistanceTO, V2)
  ) %>%
    mutate(
      NTOgr = NTO/1.15 - AirDistanceTO$Sair[1],
      NTOair = AirDistanceTO$Sair[1],
      TakeOffLength = max(BFL*as.numeric(V1 <= V2), NTO)
    )
  return(TOoutput)
  # Need to output all the dataframes!
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
    geom_hline(yintercept = inp$Srun + 200, linetype = 2) + 
    geom_line(aes(y = AccelerateStop, colour = "Accelerate-Stop")) + 
    geom_line(aes(y = AccelerateContinue, colour = "Accelerate-Continue")) +
    geom_hline(yintercept = TOoutput$BFL, linetype = 2, colour = "blue") +
    geom_vline(xintercept = TOoutput$V1, linetype = 2, colour = "blue") +
    geom_hline(yintercept = TOoutput$NTO, linetype = 2, colour = "orange") +
    geom_hline(yintercept = TOoutput$NTOgr, colour = "orange") +
    geom_vline(xintercept = TOoutput$V2, colour = "orange") +
    geom_label(aes(x = TOoutput$V2, y = 2, label = "V2")) +
    geom_label(aes(x = 10, y = TOoutput$TakeOffLength, label = "Take-off Field Length",vjust = 0.5)) +
    geom_label(aes(x = 10, y = inp$Srun, label = "Maximum Paved Runway",vjust = 0.5)) + 
    geom_label(aes(x = 10, y = inp$Srun + 200, label = "Distance to Screen",vjust = 0.5))
}

## Climb ======================================================================
ClimbFunction <- function(inp, heights) {
  #--- Add in the various interested heights
  Climboutput <- RepeatRows(inp, heights)
  Climboutput <- cbind(Climboutput, heights)
  Climboutput$type <- factor(Climboutput$type, levels = heights$type, ordered = TRUE)
  Climboutput <- StandardAtomsphere(Climboutput)
  #--- Calculate the aerodynamic parameters and Vinf values
  Climboutput <- Climboutput %>%
    mutate(
      Clmax = Clclean + flaps*Clflaps + hls*Clhls,
      Cd0 = Cd0clean + lg*Cd0lg + feathered*Cd0propfea + unfeathered*Cd0propunfea +
        (flaps + hls)*(Cd0flaps + Cdiflaps),
      Vstall = Vmin(rho, WS, Clmax),
      Vsafe = safe*Vstall,
      Vcruise = Mach*a,
      Vplot = Vcruise*1.25
    ) %>%
    rowwise() %>%
    do(data.frame(
      .,
      Vinf = c(
        seq(.$Vstall, .$Vsafe, length = 5),
        seq(.$Vsafe, .$Vcruise, length = 20),
        seq(.$Vcruise, .$Vplot, length = 20)
      ),
      Vname = c(
        "Vstall", rep("Vinf", 4),
        "Vsafe", rep("Vinf", 19),
        "Vcruise", rep("Vinf", 19))
      )) %>%
    ungroup() %>%
    data.frame(.)
  Climboutput <- Climboutput %>%
    mutate(
      qinf = 1/2 * rho * Vinf^2,
      Cl = W / (qinf * S),
      Cd = Cd0 + K * Cl^2,
      D = qinf * S * Cd,
      etaprop = etaprop(Vinf),
      PA = PA(Pshafteng, Ne, Vinf),
      TA = PA / Vinf) %>%
    rowwise() %>%
    do(data.frame(., ClimbRatesFunction(.$PA, .$Cd0, .$rho, .$Vinf, .$S, .$K, .$W))) %>%
    ungroup() %>%
    data.frame(.)
  return(Climboutput)
}

ClimbFunctionPlot <- function(inp) {
## Height values for the climb
heights <- data.frame(type = c("Takeoff", "2nd Seg", "2nd Seg OEI Feathered", "2nd Seg OEI Unfeathered", "Cruise", "Ceiling", "Landing"),
                      h = c(0, inp$Hobs, inp$Hobs, inp$Hobs, inp$AltCruise, inp$AltCeil, 0),
                      Ne = c(2, 2, 1, 1, 2, 2, 2),
                      flaps = c(1, 1, 1, 1, 0, 0, 0),
                      hls = c(0, 0, 0, 0, 0, 0, 1),
                      lg = c(1, 0, 0, 0, 0, 0, 1),
                      feathered = c(0, 0, 1, 0, 0, 0, 0),
                      unfeathered = c(0, 0, 0, 1, 0, 0, 0),
                      safe = c(1.1, 1.2, 1.2, 1.2, 1.2, 1.2, 1.3))

Climboutput <- ClimbFunction(inp, heights)

# Percentage Gradient Climb Plot
ggplot(Climboutput, aes(x=Vinf, y=PerGrad, group = type, colour = type)) +
  geom_path() +
  geom_point(aes(shape = Vname, size = ifelse(Vname == "Vinf", 0, 1))) +
  geom_hline(aes(yintercept = 1.5, colour = "2nd Seg OEI")) +
  geom_text(aes(x = min(Vinf), y = 1.5, colour = "2nd Seg OEI"),
            label = "Minimum 2nd Seg Climb OEI", hjust = 0, vjust = 1.5,
            show.legend = FALSE) +
  scale_size(range = c(0,3)) +
  scale_shape_manual(
    breaks = c("Vstall", "Vsafe", "Vcruise"),
    values = c("Vstall" = 2, "Vsafe" = 0, "Vcruise" = 1, "Vinf" = 1)) +
  guides(size = FALSE) +
  labs(list(title = "Percentage Graidents", x = "Vinf (m/s)", y = "Percentage Gradient (%)",
            colour = "Mission Segment", shape = "Velocity")) +
  ylim(-10, NA)

# Climb Rate Plot (ft/min)
ggplot(Climboutput, aes(x=Vinf, y = ClimbRate / 0.3 * 60, group = type, colour = type)) + 
  geom_path() + 
  geom_point(aes(shape = Vname, size = ifelse(Vname == "Vinf", 0, 1))) + 
  geom_hline(aes(yintercept = 100, colour = "Ceiling")) +
  geom_text(aes(x = min(Vinf), y = 100, colour = "Ceiling"),
            label = "Minimum Ceiling Rate of Climb", hjust = 0, vjust = 1.5,
            show.legend = FALSE) +
  geom_hline(aes(yintercept = 300, colour = "Cruise")) +
  geom_text(aes(x = min(Vinf), y = 300, colour = "Cruise"),
            label = "Minimum Cruise Rate of Climb", hjust = 0, vjust = 1.5,
            show.legend = FALSE) +
  scale_size(range = c(0,3)) + 
  scale_shape_manual(
    breaks = c("Vstall", "Vsafe", "Vcruise"),
    values = c("Vstall" = 2, "Vsafe" = 0, "Vcruise" = 1, "Vinf" = 1)) +
  guides(size = FALSE) +
  labs(list(title = "Climb Rates (Vv)", x = "Vinf (m/s)", y = "Climb Rate (ft/min)", 
            colour = "Mission Segment", shape = "Velocity")) +
  ylim(-10, NA)

# n Loading
ggplot(Climboutput, aes(x=Vinf, y=nload, group = type, colour = type)) +
  geom_path() +
  geom_point(aes(shape = Vname, size = ifelse(Vname == "Vinf", 0, 1))) +
  scale_size(range = c(0,3)) + 
  scale_shape_manual(
    breaks = c("Vstall", "Vsafe", "Vcruise"),
    values = c("Vstall" = 2, "Vsafe" = 0, "Vcruise" = 1, "Vinf" = 1)) +
  guides(size = FALSE) +
  labs(list(title = "n Loading", x = "Vinf (m/s)", y = "n Loading",
            colour = "Mission Segment", shape = "Velocity"))

}


## Landing ======================================================================
LandingLength <- function(inp, V2 = 1.1 * inp$VsLD) {
  #--- Determine the distance required for landing
  AirDistLD <- inp
  AirDistLD$type <- "All Engines"
  AirDistLD$Ne <- 2
  AirDistLD <- AirDistLD %>%
    mutate(h = Hobsland) %>%
    StandardAtomsphere(.) %>%
    mutate(
      Clapp = Clclean + Clhls,
      Cd0 = Cd0clean + Cd0lg + Cd0flaps + Cdiflaps,
      Vapp = 1.3* VsLD,
      qinf = 1/2 * rho * Vapp^2,
      etaprop = etaprop(Vapp),
      PA = PA(Pshafteng, Ne, Vapp) * 0.05,
      gammadeg = max(-3, ClimbRatesFunction(PA, Cd0, rho, Vapp, S, K, W)[[1]]),
      gamma = gammadeg * pi/180,
      L = W * cos(gamma),
      Cl = L / (qinf * S),
      Cd = Cd0 + K * Cl^2,
      D = qinf * S * Cd,
      TR = D - W * sin(gamma),
      PR = TR * Vapp / etaprop,
      R = Vapp^2 / (0.2 * g),
      SF = R * sin(gamma),
      hF = R * (1 - cos(gamma)),
      SA = (Hobsland - hF) / tan(gamma),
      Sair = ifelse(SA >0, SA + SF, sqrt(R^2 - (R-Hobsland)^2)))
  #--- Set up the initial parameters to solve for
  LD <- inp
  LD$segment <- "Landing"
  LD$type <- "Landing"
  LD$Ne <-0
  LD$mu <- c(
    as.double(groundmu["Dry Concrete", "brakeson"])
  )
  #---- Determine aerodynamic parameters
  LD <- mutate(LD, h = 0) %>%
    StandardAtomsphere(.) %>%
    mutate(
      K = Keff(K, hground, b),
      Cd0 = Cd0clean + Cd0lg + Cd0flaps + Cdiflaps,
      Cd0 = Cd0 + as.numeric(Ne == 1) * Cd0propfea + as.numeric(Ne == 0) * 2 * Cd0propunfea,
      Cl = Cl0 + Clflaps,
      Cd = Cd0 + K * Cl^2
    )
  #--- Output the results
  LDoutput <- data.frame(
    Vapp = AirDistLD$Vapp,
    VTD = 1.15*inp$VsLD
  ) %>%
    mutate(
      LD = (VTD * 3 + 
        integrate(function(V) V/GroundAcceleration(LD, V, distancecalc = TRUE),
                   lower = 1.15*inp$VsLD, upper = 0)[[1]] +
          AirDistLD$Sair[1]) * 1.67,
      LDgr = LD/1.67 - AirDistLD$Sair[1],
      LDair = AirDistLD$Sair[1]
    )
  return(LDoutput)
  # Need to output all the dataframes!
}

