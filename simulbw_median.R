phi.hat.11.med = function(x, p, T){
  nT = floor(T^0.66)
  W = 1 # any diagonal matrix
  res = 0
  for(r in 1:p){
    res.in = array(NA, dim = floor(T/nT))
    for(j in 1:floor(T/nT)){
      term1 = (sigma.hat(x, r, (j*nT + 1)/T, T))^2
      term2 = (1 - a.hat(x, r, (j*nT + 1)/T, T))^(-2)
      res.in[j] = term1*term2
      # cat("In phi.hat.11:\n")
      # cat("term1 =", term1, "\nterm2 =", term2, "\nres.in =", res.in, "\n")
    }
    res = res + W[r][r] * ((median(res.in))^(-2))
  }
  res*(delta.barhat.121(T)^2)*((2*pi)^2)
}

phi.hat.12.med = function(x, p, T){
  nT = floor(T^0.66)
  W = 1 # same as phi.hat.11
  res = 0
  for(r in 1:p){
    up = array(NA, dim = floor(T/nT))
    down = array(NA, dim = floor(T/nT))
    for(j in 1:floor(T/nT)){
      sigma = sigma.hat(x, r, (j*nT + 1)/T, T)
      a = a.hat(x, r, (j*nT + 1)/T, T)
      #cat("(j*nT + 1)/T =", (j*nT + 1)/T, "\n")
      #cat("a =", a, "sigma =", sigma, "\n")
      #cat("true a =", 0.8*(cos(1.5*pi/180)+cos(4*pi*(j*nT+1)/T)), "true sigma = 1\n")
      up[j] = (sigma^2)*a*((1 - a)^(-4))
      down[j] = (sigma^2)*((1 - a)^(-2))
    }
    res = res + W[r][r]*((median(up)/median(down))^2)
  }
  
  res*36
}

# merge everything together
bw_opt_med = function(x){
  T = length(x)
  p = 1
  
  # optimal bandwidth
  phi11 = phi.hat.11.med(x, p, T)
  phi12 = phi.hat.12.med(x, p, T)
  #cat("phi11 =", phi11, "phi12", phi12, "\n")
  phi1 = phi11/(phi12^5)
  phi2 = phi12/(phi11^5)
  b1 = 0.46*(phi1^(1/24))*T^(-1/6)
  b2 = 3.56*(phi2^(1/24))*T^(-1/6)
  #cat("phi11 =", phi11, "phi12 =", phi12, "\n")
  #cat("phi1 =", phi1, "phi2 =", phi2, "\n")
  #cat("b1 =", b1, "b2 =", b2, "\n")
  c(Re(b1), Re(b2))
}

# DK-HAC
JHat1MedBW = function(x, b1=NA, b2=NA){
  if(is.na(b1) | is.na(b2)){
    b = bw_opt_med(x)
    b1 = b[1]
    b2 = b[2]
  }
  T = length(x)
  nT = T^(0.66)
  b = floor(1/b1)
  if(b > T)
    b = T-1
  gamma = rep(0, b+1)
  
  for(k in 0:b){
    for(r in 0:floor((T-nT)/nT)){
      s = (k+1):T
      cT = sum(x[s]*x[s-k]*k2(((r+1)*nT - (s-k/2))/(T*b2)))/(T*b2)
      gamma[k+1] = gamma[k+1] + cT
    }
  }
  gamma = gamma*nT/(T-nT)
  #cat("gamma =", gamma[1:5], "\n")
  ker = k1(b1*(0:b))
  res = gamma*ker
  res[1]+2*sum(res[-1])
}

JHat1MedcT = function(x, b1=NA, b2=NA){
  if(is.na(b1) | is.na(b2)){
    b = bw_opt(x)
    b1 = b[1]
    b2 = b[2]
  }
  T = length(x)
  nT = T^(0.66)
  b = floor(1/b1)
  if(b > T)
    b = T-1
  gamma = rep(0, b+1)
  
  for(k in 0:b){
    cT = array(NA, dim = floor((T-nT)/nT + 1))
    for(r in 0:floor((T-nT)/nT)){
      s = (k+1):T
      cT[r+1] = sum(x[s]*x[s-k]*k2(((r+1)*nT - (s-k/2))/(T*b2)))/(T*b2)
    }
    gamma[k+1] = median(cT)
  }
  #cat("gamma =", gamma[1:5], "\n")
  ker = k1(b1*(0:b))
  res = gamma*ker
  res[1]+2*sum(res[-1])
}

JHat1Med = function(x, b1=NA, b2=NA){
  if(is.na(b1) | is.na(b2)){
    b = bw_opt_med(x)
    b1 = b[1]
    b2 = b[2]
  }
  T = length(x)
  nT = T^(0.66)
  b = floor(1/b1)
  if(b > T)
    b = T-1
  gamma = rep(0, b+1)
  
  for(k in 0:b){
    cT = array(NA, dim = floor((T-nT)/nT + 1))
    for(r in 0:floor((T-nT)/nT)){
      s = (k+1):T
      arr = x[s]*x[s-k]
      #cT[r+1] = median(arr) * length(arr)/(T*b2)
      gamma[k+1] = gamma[k+1] + (T^(1/6))*sum(median(arr)*k2(((r+1)*nT - (s-k/2))/(T*b2)))/(T*b2)
    }
    #cat("cT =", cT, "\n")
    #gamma[k+1] = median(cT)
  }
  gamma = gamma*nT/(T-nT)
  #cat("gamma =", gamma[1:5], "\n")
  #plot(1:length(gamma), gamma)
  ker = k1(b1*(0:b))
  res = gamma*ker
  res[1]+2*sum(res[-1])
}

