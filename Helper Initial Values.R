#----------------------------
#--- Initial Values
#============================

## Initial Inputs ======================================================================
input_initial <- data.frame(
  # Configuration
  S = 24.5,
  b = 22.1359,
  AR = 20,
  Taper = 0.4,
  
  # Aircraft Efficiency
  e = 0.79929,
  K = 0.01991,
  
  # Aerodynamic Properties
  Clclean = 1.8,
  Clflaps = 1.1,
  Clhls = 1.1, 
  ClG = 0.31801,
  Cd0clean = 0.02291,
  Cd0lg = 0.02857,
  Cd0OEI = 0.00580,
  Cd0flaps = 0.03885,
  Cdiflaps = 0.09486,
  Cd0propfea = 0.00580,
  Cd0propunfea = 0.04636,
  hground = 2.5,
  
  # Weights
  m = 6900,
  W = 67664.85,
  WS =  2761.831,
  
  # Shaft Power per Engin
  Pshafteng = 668100,
  Pshaft = 1336200,
  etamech = 0.91238
)

## Specifications ======================================================================
specs_decript <- data.frame(
  Variable = c("Mp", "Wp", "E", "Dens", "Srun", "Hobs", "Vappmax", "PerGrad2Seg", "ClimbCruise", "ClimbCeil", "AltCruise", "AltCeil", 
               "Mach", "Range", "LoadMax", "LoadMin"),
  Description = c("Payload Mass (kg)","Payload Weight (N)","Battery Specific Energy (J/kg)","Battery Density (kg/m^3)",
                  "Runway Length (m)","Screen Height (m)","Maximum Landing Speed (m/s)","Climb Gradient 2nd Segment (%)",
                  "Climb at Cruise (m/s)","Climb at Ceiling (m/s)","Altitude of Cruise (m)","Altitude of Ceiling (m)",
                  "Cruise Mach Number (M)","Range (m)","Load Factor Max (n)","Load Factor Min (n)"),
  Value = c(120 * 6, 120 * 6 * 9.8065, 1e6, 2750, 1000, 35 * 0.3048, 100 * 0.5144, 1.5, 300 * 0.3048 / 60,
            100 * 0.3048 / 60, 10000 * 0.3048, 12000 * 0.3048, 0.25, 1000e3, 3.5, -1.5)
)

specs <- setNames(data.frame(t(specs_decript["Value"])),  t(specs_decript["Variable"]))

## Desnities ======================================================================
specs$rhoCruise <- StandardAtomsphere(data.frame(h = specs$AltCruise))$rho
specs$rhoCeil <- StandardAtomsphere(data.frame(h = specs$AltCeil))$rho

## inp Dataframe ======================================================================
inp = data.frame(segment = factor(
  "NA",
  levels = c(
    "Takeoff",
    "Transition",
    "Segment 1",
    "Segment 2",
    "Segment 3",
    "Segment 4",
    "Cruise",
    "Descend",
    "Flare",
    "Landing",
    "NA"
  ),
  ordered = TRUE
),
input_initial,
specs)
rownames(inp) <- NULL
