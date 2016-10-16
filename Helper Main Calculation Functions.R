#----------------------------
#--- Main Calculation Functions
#============================

## Groundroll ======================================================================
GroundAcceleration <- function(inp) {
  
}

## Takeoff ======================================================================
# TakeOff <- function(inp) {
  # Set up the initial parameters to solve for
  TO <- RepeatRows(inp, 3)
  TO$segment <- "Takeoff"
  TO$type <- c("All Engines", "One Engine Down", "Rejected Take-Off")
  TO$Ne <- c(2, 1, 0)
  TO$mu <- c(
    as.double(groundmu["Dry Concrete", "brakesoff"]),
    as.double(groundmu["Dry Concrete", "brakesoff"]),
    as.double(groundmu["Dry Concrete", "brakeson"])
  )
  # Determine aerodynamic parameters
  TO <- mutate(TO, h = 0) %>%
    StandardAtomsphere(.) %>%
    mutate(
      Cd0 = Cd0clean + Cd0lg + Cd0flaps + Cdiflaps,
      Cd0 = Cd0 + (2 - Ne) * Cd0propunfea,
      Cl = Cl0 + Clflaps
    )
# }