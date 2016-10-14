#----------------------------
#--- Numerical Methods
#============================

trapzd <- function(func, a, b, s, n) {
  if (n == 1) {
    s = 1/2 * (b - a) * (func(a) + func(b))
    ssum = s
  }
  else {
    it = 2^(n-2)
    h = (b - a) / it
    x = a + 0.5 * h
    ssum = 0
    for (j in 1:it) {
      ssum = ssum + func(x)
      x = x + h
    }
    ssum = 0.5 * (s + h*ssum)
  }
  return(ssum)
}

RecurTrapzd <- function(func, a, b, toler) {
  j = 0
  s = 0
  s_new = 1.2345
  while (abs(s_new - s) > toler) {
    j = j + 1
    s = s_new
    s_new = trapzd(func, a, b, s, j)
  }
  return(s_new)
}

ModifiedSecant <- function(func, xr, del, toler, positive = FALSE) {
  xr_old = xr - 1
  while (abs(xr - xr_old) > toler) {
    fx1 = func(xr)
    fx2 = func(xr + del)
    xr_old = xr
    xr = xr_old - (del * fx1) / (fx2 - fx1)
    if (xr < 0 & positive == TRUE)
      return(NA)
  }
  return(xr)
}

fx <- function(x) x^2 - x - 7

