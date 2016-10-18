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