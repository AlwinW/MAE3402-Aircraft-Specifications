#----------------------------
#--- Energy Usage
#============================

## Energy Summary ======================================================================
EnergyUsage <- function(TO, AirDistanceTO, LD, AirDistanceLD, inp, reserve = 0.05, sf = 1.05) {
## Take Off ======================================================================
  #--- Clean up input
  TOall <- filter(TO, type == "All Engines")
  AirDistanceall <- filter(AirDistanceTO, type == "All Engines")
  #--- Take-off Ground Roll
  TOgr <- integrate(function(V)
      GroundAcceleration(TOall, V, energycalc = TRUE),
    lower = 0, upper = 1.1*inp$VsTO
  )[[1]] +
    PA(inp$Pshafteng, 2, 1.1 * inp$VsTO) * 3
  #--- Take-off Transition
  TOair <- integrate(Vectorize(function(h)
    PA(inp$Pshafteng, 2, 1.15*inp$VsTO) / 
      (sin(AirDistanceall$gamma) * 1.15*inp$VsTO)),
    lower = 0, upper = AirDistanceall$hTR
  )[[1]]
## Climb Segments ======================================================================
  #--- 1st Semgnet Climb
  if (AirDistanceall$SC > 0) {
    Seg1 <- integrate(Vectorize(function(h)
      PA(inp$Pshafteng, 2, 1.15*inp$VsTO) / 
        (sin(atan(AirDistanceall$hTR/(AirDistanceall$Sair - AirDistanceall$ST))) * 1.15*inp$VsTO)),
      lower = AirDistanceall$hTR, upper = inp$Hobs
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
    lower =  max(AirDistanceall$hTR, inp$Hobs), upper = inp$AltFlaps
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
  #--- Time Required for Desc1 from cruise to 400ft
  Desc1inp <- mutate(inp, Cd0 = Cd0clean + 2*Cd0propfea)
  Desc1time = integrate(function(h)
    DescentZeroAccel(Desc1inp, h, inp$Vcruise, timecalc = TRUE),
    lower = inp$AltCruise, upper = inp$AltFlaps
  )[[1]]
  #--- Approximate Power Required (5% Reserve)
  Desc1 = (Desc1time) * PA(inp$Pshafteng, 2, inp$Vcruise)/etaprop(inp$Vcruise) * reserve
  #--- Flaps
  Desc2inp <- inp %>%
    mutate(
      h = AltFlaps,
      Cd0 = Cd0clean + Cd0flaps + Cdiflaps,
      Ne = 0
    )  %>%
    StandardAtomsphere(.)
  Desc2time <- integrate(function(V)
    1/ AccelerateEnergy(Desc2inp, V, distancecalc = TRUE),
    lower =  1.3*inp$VsLD, upper = inp$Mach*Desc2inp$a)[[1]]
  #--- Approximate Power Required (5% Reserve)
  Desc2 = (Desc2time) * PA(inp$Pshafteng, 2, 1.3*inp$VsLD)/etaprop(1.3*inp$VsLD) * reserve
  #--- Time Required for Desc3 from 400ft to 50ft
  Desc3inp <- mutate(inp, Cd0 = Cd0clean + 2*Cd0propfea + Cd0flaps + Cdiflaps)
  Desc3time = integrate(function(h)
    DescentZeroAccel(Desc3inp, h, 1.3*inp$VsLD, timecalc = TRUE),
    lower = inp$AltFlaps, upper = inp$Hobsland
  )[[1]]
  #--- Approximate Power Required (5% Reserve)
  Desc3 = (Desc3time) * PA(inp$Pshafteng, 2, 1.3*inp$VsLD)/etaprop(inp$Vcruise) * reserve
  
## Landing ======================================================================
  #--- Take-off Ground Roll
  LDgr <- integrate(function(V)
    GroundAcceleration(LD, V, energycalc = TRUE),
    lower = 0, upper = 1.15*inp$VsLD)[[1]] +
    PA(inp$Pshafteng, 0, 1.15 * inp$VsLD) * 3
  #--- Take-off Transition
  LDair <- integrate(Vectorize(function(h)
    AirDistanceLD$PR/ 
      (sin(AirDistanceLD$gamma) * 1.15*inp$VsLD)),
    lower = AirDistanceLD$Hobsland, upper = 0)[[1]]
  
## Distances ======================================================================
  Distances <- data.frame(
    TOgr = integrate(function(V) V/GroundAcceleration(filter(TO, type == "All Engines"), V, distancecalc = TRUE), 
                     lower = 0, upper = 1.1*inp$VsTO)[[1]] + 3 * 1.1*inp$VsTO,
    TOair = AirDistanceall$Sair,
    
    Seg1 = 0,
    Seg2 = integrate(function(h) ClimbEnergy(Seg2inp, h, distancecalc = TRUE),
      lower =  max(AirDistanceall$hTR, inp$Hobs), upper = inp$AltFlaps)[[1]],
    Seg3 = integrate(function(V) V/AccelerateEnergy(Seg3inp, V, distancecalc = TRUE), 
                     lower = 1.2*inp$VsTO, upper = inp$Mach*Seg3inp$a)[[1]],
    Seg4 = integrate(function(h) ClimbEnergy(Seg4inp, h, distancecalc = TRUE),
                     inp$AltFlaps, upper = inp$AltCruise)[[1]],
    
    Cruise = 0,
    
    Desc1 = integrate(function(h) DescentZeroAccel(Desc1inp, h, inp$Vcruise, distancecalc = TRUE),
                     lower = inp$AltCruise, upper = inp$AltFlaps)[[1]],
    Desc2 = integrate(function(V) V/AccelerateEnergy(Desc2inp, V, distancecalc = TRUE),
                     lower =  1.3*inp$VsLD, upper = inp$Mach*Desc2inp$a)[[1]],
    Desc3 = integrate(function(h) DescentZeroAccel(Desc3inp, h, 1.3*inp$VsLD, distancecalc = TRUE),
                      lower = inp$AltFlaps, upper = inp$Hobsland)[[1]],
    
    LDair = AirDistanceLD$Sair,
    LDgr = 1.15*inp$VsLD * 3 + 
                      integrate(function(V) V/GroundAcceleration(LD, V, distancecalc = TRUE),
                      lower = 1.15*inp$VsLD, upper = 0)[[1]] +
                      AirDistanceLD$Sair[1],
    
    Systems = 0
  )
  Distances$Cruise = inp$Range - sum(Distances)
  Distances$Systems = inp$Range
  
## Cruise ======================================================================
  Cruiseinp <- mutate(inp, h = AltCruise) %>%
    StandardAtomsphere(.) %>%
    mutate(
      Vinf = Vcruise,
      qinf = 1/2 * rho * Vinf^2,
      Cl = W/(qinf * S),
      Cd0 = Cd0clean,
      Cd = Cd0 + K * Cl^2,
      D = qinf * S * Cd,
      TR = D,
      PR = TR * Vinf,
      etaprop = etaprop(Vinf)
      )
  Cruise = Cruiseinp$TR / Cruiseinp$etaprop * Distances$Cruise
  Cruiseonly = Cruiseinp$TR / Cruiseinp$etaprop * inp$Range

## Systems ======================================================================
  ## TO DO ####
  Systems = 0
  
  EnergySummary <- data.frame(
    TOgr = TOgr,
    TOair = TOair,
    Seg1 = Seg1,
    Seg2 = Seg2,
    Seg3 = Seg3,
    Seg4 = Seg4,
    Cruise = Cruise,
    Desc1 = Desc1,
    Desc2 = Desc2,
    Desc3 = Desc3,
    LDair = LDair,
    LDgr = LDgr,
    Systems = Systems
  )
  
  Mbatteries = sum(energysummary/inp$E/inp$etamech*sf)
  MbatteriesEst = sum(Cruiseonly/inp$E/inp$etamech*sf)
  
  return(list(
    EnergySummary= EnergySummary,
    Mbatteries = Mbatteries,
    MbatteriesEst = MbatteriesEst,
    Distances = Distances
    ))
}

# EnergySummaryAll <- EnergyUsage(TO, AirDistanceTO, LD, AirDistanceLD, inp, reserve = 0.1, sf = 1.05)

EnergyUsagePlot <- function(EnergySummaryAll) {
  #--- Energy
  EnergySummaryPlotData = EnergySummaryAll$EnergySummary %>%
    gather(key = Segment, value = EnergyUsage)
  EnergySummaryPlotData$Segment <- factor(
    EnergySummaryPlotData$Segment, levels = EnergySummaryPlotData$Segment, ordered = TRUE)
  
  ggplot(data = EnergySummaryPlotData, aes(x = Segment)) + 
    geom_point(aes(y = EnergyUsage))
  
  #--- Average Power
  AveragePowerPlotData = EnergySummaryAll$EnergySummary/EnergySummaryAll$Distances 
  AveragePowerPlotData <- gather(AveragePowerPlotData, key = Segment, value = AveragePower)
  AveragePowerPlotData$Segment <- factor(
    AveragePowerPlotData$Segment, levels = AveragePowerPlotData$Segment, ordered = TRUE)
  
  ggplot(data = AveragePowerPlotData, aes(x = Segment)) + 
    geom_point(aes(y = AveragePower))
  
  #--- Battery Breakdown
  BatteryFracPlotData = EnergySummaryAll$EnergySummary / as.numeric(EnergySummaryAll$Mbatteries)
  BatteryFracPlotData <- gather(BatteryFracPlotData, key = Segment, value = BatteryFraction)
  BatteryFracPlotData$Segment <- factor(
    BatteryFracPlotData$Segment, levels = BatteryFracPlotData$Segment, ordered = TRUE)
  
  ggplot(data = BatteryFracPlotData, aes(x = Segment)) + 
    geom_point(aes(y = BatteryFraction))
}






