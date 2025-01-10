k1 = function(x){
  if(x == 0)
    return(1)
  res = (25/(12*(pi^2)*(x^2))) * (sin(6*pi*x/5)/(6*pi*x/5) - cos(6*pi*x/5))
  return(res)
}

k2 = function(x){
  if(x < 0 | x > 1)
    return(0)
  res = 6*x*(1 - x)
  return(res)
}

dtHat = function(u ,k){
  S.omega = c(-pi, -3, -2, -1, 0, 1, 2, 3, pi)
  res = 0
  i = complex(imaginary = 1)
  for(s in S.omega){
    n1 = (3/pi)*((1 + 0.8*(cos(1.5) + cos(4*pi*u))*exp(-i*s))^(-4))
    n2 = 0.8*(-4*pi*sin(4*pi*u))*exp(-i*s)
    n3 = (-1/pi)*((abs(1 + 0.8*(cos(1.5) + cos(4*pi*u))*exp(-i*s)))^(-3))
    n4 = 0.8*(-16)*(pi^2)*cos(4*pi*u)*exp(-i*s)
    res = res + (n1*n2 + n3*n4)*exp(-i*k*s)
  }
  return(res/length(S.omega))
}

dtBarHat = function(nT, T){
  res = 0
  for(k in -floor(T^(-1/6)):floor(T^(1/6))){
    res.k = 0
    for(j in 0:floor(T/nT)){
      res.k = res.k + dtHat(j*nT/T, k)
    }
    res = res + res.k*nT/T
  }
  return(res)
}

a1Hat = function(x, u, nT, T){
  t = floor(u*T)
  up = 0
  down = 0
  for(j in (t-nT+1):t){
    up = up + x[j]*x[j-1]
    down = down + (x[j-1])^2
  }
  return(up/down)
}

sigmaHat = function(x, u, nT, T){
  res = 0
  t = floor(u*T)
  for(j in (t-nT+1):t){
    res = res + (x[j] - a1Hat(x, u, nT, T) * x[j-1])^2
  }
  return(res^(1/2))
}

phiHat11 = function(x, nT, T){
  up = ((4*pi)^(-2)) * (dtBarHat(nT, T))^2
  down = 0
  for(j in 1:floor(T/nT - 1)){
    sigma = sigmaHat(x, (j*nT+1)/T, nT, T)
    a = a1Hat(x, (j*nT+1)/T, nT, T)
    down = down + (sigma^2)*((1 - a)^(-2))
  }
  down = down*nT/T
  down = down^2
  return(up/down)
}

phiHat12 = function(x, nT, T){
  up = 0
  down = 0
  for(j in 1:floor(T/nT - 1)){
    sigma = sigmaHat(x, (j*nT+1)/T, nT, T)
    a = a1Hat(x, (j*nT+1)/T, nT, T)
    up = up + ((sigma^2)*a)*((1 - a)^(-4))
    down = down + (sigma^2)*((1 - a)^(-2))
  }
  return(36*((up/down)^(2)))
}

BWselection = function(x, nT, T){
  phi11 = phiHat11(x, nT, T)
  phi12 = phiHat12(x, nT, T)
  cat("phi11 =", phi11, "phi12", phi12, "\n")
  phi1 = phi11/(phi12^5)
  phi2 = phi12/(phi11^5)
  b1 = 0.46*(phi1^(1/24))*(T^(-1/6))
  b2 = 3.56*(phi2^(1/24))*(T^(-1/6))
  return(c(Re(b1), Re(b2)))
}

JHat = function(x, nT, T, b1 = NA, b2 = NA){
  if(is.na(b1) | is.na(b2)){
    b = BWselection(x, nT, T)
    b1 = b[1]
    b2 = b[2]
    cat("b1 =", b1, "b2 =", b2, "\n")
  }
  res = 0
  for(k in (-T+1):(T-1)){
    Gamma = 0
    for(r in 0:floor((T-nT)/nT)){
      cT = 0
      if(k >= 0){
        sum = 0
        for(s in (k+1):T){
          sum = sum + x[s]*x[s-k]*k2(((r+1)*nT - (s-k/2))/(T*b2))
        }
        cT = cT + sum/(T*b2)
      }
      else{
        sum = 0
        for(s in (-k+1):T){
          sum = sum + x[s+k]*x[s]*k2(((r+1)*nT - (s+k/2))/(T*b2))
        }
        cT = cT + sum/(T*b2)
      }
      Gamma = Gamma + cT
    }
    Gamma = k1(b1*k)*Gamma*nT/(T-nT)
    res = res + Gamma
  }
  return(res)
}



