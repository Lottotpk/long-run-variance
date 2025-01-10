k1 = function(x){ # QS
  idx = which(x == 0)
  if(length(idx) != 0)
    idx2 = (1:length(x))[-idx]
  else
    idx2 = 1:length(x)
  x[idx] = 1
  x[idx2] = (25/(12*(pi^2)*(x[idx2]^2))) * (sin(6*pi*x[idx2]/5)/(6*pi*x[idx2]/5) - cos(6*pi*x[idx2]/5))
  x
}

k2 = function(x){
  idx = which(x >= 0 & x <= 1)
  if(length(idx) != 0)
    idx2 = (1:length(x))[-idx]
  else
    idx2 = 1:length(x)
  x[idx2] = 0
  x[idx] = 6*x[idx]*(1-x[idx])
  x
}

a.hat = function(x, r, u, T){
  nT = floor(T^0.66)
  t = u*T
  start = t-nT+1
  up = 0
  down = 0
  
  if(start < 1)
    return(0)
  for(j in abs(start):t){
    if(j == 1) next
    up = up + x[j][r]*x[j-1][r]
    down = down + (x[j-1][r])^2
  }
  up/down
}

sigma.hat = function(x, r, u, T){
  #cat("In sigma.hat:\n")
  nT = floor(T^0.66)
  t = u*T
  start = t-nT+1
  res = 0
  
  if(start < 1)
    return(0)
  for(j in abs(start):t){
    if(j == 1) next
    res = res + (x[j][r] - a.hat(x, r, u, T)*x[j-1][r])^2
  } 
  #cat("res =", res, "\n")
  res^(1/2)
}

phi.hat2 = function(x, p, T){
  nT = T^(0.66)
  W = 1
  res = 0
  for(r in 1:p){
    up = 0
    down = 0
    for(j in 0:floor(T/nT-1)){
      sigma = sigma.hat(x, r, (j*nT+1)/T, T)
      a = a.hat(x, r, (j*nT+1)/T, T)
      num_up = (sigma*a)^2
      den_up = (1-a)^4
      up = up + num_up/den_up
      
      num_down = sigma^2
      den_down = (1-a)^2
      down = down + num_down/den_down
    }
    up = 18 * W[r][r] * ((up*nT/T)^2)
    down = W[r][r] * ((down*nT/T)^2)
    
    res = res + up/down
  }
  res
}

D1.hat = function(u, k){
  S.omega = c(-pi, -3, -2, -1, 0, 1, 2, 3, pi)
  n = length(S.omega)
  i = complex(imaginary = 1)
  sum = 0
  
  for(n.omega in 1:n){
    e = exp(-1*i*S.omega[n.omega])
    term1 = (3/pi) * (1 + 0.8*(cos(1.5) + cos(4*pi*u))*e)^(-4)
    term2 = 0.8*(-4)*pi*sin(4*pi*u)*e
    term3 = (1/pi) * (abs(1 + 0.8*(cos(1.5) + cos(4*pi*u))*e))^(-3)
    term4 = 0.8*(-16)*(pi^2)*cos(4*pi*u)*e
    sum = sum + (term1*term2 - term3*term4)*exp(i*k*S.omega[n.omega])
  }
  (sum/n)^2
}

cT.hat = function(x, u, k, T){
  
}

D2.hat = function(u, p, T){
  range = floor(T^(4/25))
  res = 0
  for(r in 1:p){
    for(l in -range:range){
      
    }
  }
}

Jseq = function(x, b1=NA, b2=NA){
  if(is.na(b1) | is.na(b2)){
    b = bw_opt(x)
    b1 = b[1]
    b2 = b[2]
  }
  #cat("JHat1 :", "b1 =", b1, "b2 =", b2, "\n")
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