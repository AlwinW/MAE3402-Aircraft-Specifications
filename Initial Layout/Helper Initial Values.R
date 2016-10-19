#----------------------------
#--- Initial Values
#============================

## Initial Inputs ======================================================================
input_initial <- data.frame(
  S = 24.5,
  b = 22.1359,
  AR = 20,
  e = 0.79929,
  K = 0.01991,
  
  Cd0 = 0.02291,
  Clclean = 1.8,
  Clflaps = 0.9169118,
  Clhls = 0.9169118,
  
  m = 6800,
  W = 66684.2,
  WS = 2721.804,
  
  P0eng = 705000,
  P0 = 1410000,
  
  ClG = 0.25,
  Cd0G = 0.034365,
  hground = 2.5
)

resolution = 10
inputvals <- input_initial

## Specifications ======================================================================
specifications_description <- data.frame(
  Variable = c("Mp", "Wp", "E", "Dens", "Srun", "Hobs", "Vappmax", "PerGrad2Seg", "ClimbCruise", "ClimbCeil", "AltCruise", "AltCeil", 
               "Mach", "Range", "LoadMax", "LoadMin"),
  Description = c("Payload Mass (kg)","Payload Weight (N)","Battery Specific Energy (J/kg)","Battery Density (kg/m^3)",
                  "Runway Length (m)","Screen Height (m)","Maximum Landing Speed (m/s)","Climb Gradient 2nd Segment (%)",
                  "Climb at Cruise (m/s)","Climb at Ceiling (m/s)","Altitude of Cruise (m)","Altitude of Ceiling (m)",
                  "Cruise Mach Number (M)","Range (m)","Load Factor Max (n)","Load Factor Min (n)"),
  Value = c(120 * 6, 120 * 6 * 9.8065, 1e6, 2750, 1000, 35 * 0.3048, 100 * 0.5144, 1.5, 300 * 0.3048 / 60,
            100 * 0.3048 / 60, 10000 * 0.3048, 12000 * 0.3048, 0.25, 1000e3, 3.5, -1.5)
)

specifications  <- t(specifications_description["Value"])
colnames(specifications) <- t(specifications_description["Variable"])

