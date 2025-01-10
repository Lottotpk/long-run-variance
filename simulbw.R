# AR(1)
x = arima.sim(model=list(ar=.5), n=1000)

# Classical lrv estimation
JCla = function(x, b=NA){
  if(is.na(b)){
    b = ceiling(2*length(x)^(1/3))
  }
  cHat = c(acf(x, type="cov", lag.max=b, plot=FALSE)$acf)	
  kernel = 1-(0:b)/b
  temp = cHat*kernel
  temp[1] + 2*sum(temp[-1])
}


## This article estimation
# Optimal kernels
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

# Estimation functions
a.hat = function(x, r, u, T){
  nT = floor(T^0.66)
  t = u*T
  start = t-nT+1
  up = 0
  down = 0
  
  if(start < 1)
    return(0)
  for(j in start:t){
    if(j == 1) next
    if(j >= T) next
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
    return(x[t]^(1/2))
  for(j in start:t){
    if(j == 1) next
    if(j >= T) next
    res = res + (x[j][r] - a.hat(x, r, u, T)*x[j-1][r])^2
  } 
  #cat("res =", res, "\n")
  res^(1/2)
}       

delta.hat.121 = function(u, k){
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
  sum/n
}

delta.barhat.121 = function(T){
  nT = floor(T^0.66)
  res = 0
  T_floor = floor(T^(-1/6))
  # cat("In delta.barhat.121:\n")
  # cat("nT =", nT, "T_floor", T_floor, "\n")
  for(k in -T_floor:T_floor){
    res.in = 0
    for(j in 0:floor(T/nT)){
      res.in = res.in + delta.hat.121(j*nT/T, k)
      # cat("res.in =", res.in, "\n")
    }
    res = res + ((nT/T) * res.in)
  }
  res
}

phi.hat.11 = function(x, p, T){
  nT = floor(T^0.66)
  W = 1 # any diagonal matrix
  res = 0
  for(r in 1:p){
    res.in = 0
    for(j in 1:floor(T/nT)){
      term1 = (sigma.hat(x, r, (j*nT + 1)/T, T))^2
      term2 = (1 - a.hat(x, r, (j*nT + 1)/T, T))^(-2)
      res.in = res.in + term1*term2
      # cat("In phi.hat.11:\n")
      # cat("term1 =", term1, "\nterm2 =", term2, "\nres.in =", res.in, "\n")
    }
    res = res + W[r][r] * ((nT*res.in/T)^(-2))
  }
  res*(delta.barhat.121(T)^2)*((2*pi)^2)
}

phi.hat.12 = function(x, p, T){
  nT = floor(T^0.66)
  W = 1 # same as phi.hat.11
  res = 0
  for(r in 1:p){
    up = 0
    down = 0
    for(j in 1:floor(T/nT)){
      sigma = sigma.hat(x, r, (j*nT + 1)/T, T)
      a = a.hat(x, r, (j*nT + 1)/T, T)
      #cat("(j*nT + 1)/T =", (j*nT + 1)/T, "\n")
      #cat("a =", a, "sigma =", sigma, "\n")
      #cat("true a =", 0.8*(cos(1.5*pi/180)+cos(4*pi*(j*nT+1)/T)), "true sigma = 1\n")
      up = up + (sigma^2)*a*((1 - a)^(-4))
      down = down + (sigma^2)*((1 - a)^(-2))
    }
    res = res + W[r][r]*((up/down)^2)
  }
  
  res*36
}

# merge everything together
bw_opt = function(x){
  T = length(x)
  p = 1
  
  # optimal bandwidth
  phi11 = phi.hat.11(x, p, T)
  phi12 = phi.hat.12(x, p, T)
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
JHat1 = function(x, b1=NA, b2=NA){
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

# try
x = arima.sim(model = list(ar = .5), n = 1000)
T = length(x)
nT = floor(T^(0.66))
#c(BWselection(x, nT, T), bw_opt(x), 0.25*T^(-1/6), 0.45*T^(-1/6))
c(JHat1(x), JHat1(x, b1=0.25*T^(-1/6), b2=0.45*T^(-1/6)), JCla(x))

# timer
start = Sys.time()
n.all = seq(100, 1000, 100)
res1 = array(NA, length(n.all))
res2 = array(NA, length(n.all))
res3 = array(NA, length(n.all))
for(i in n.all){
  x=arima.sim(model=list(ar=0.5), n=i)
  res1[i/100] = JHat1Med(x)
  res2[i/100] = JHat1(x)
  res3[i/100] = JCla(x)
}
plot(n.all, res1, col="red", ylim = c(0, 10))
points(n.all, res2, col="blue")
points(n.all, res3, col="green")
abline(h = 4)
#x[sample(length(x), ceiling(0.5*length(x)^(1/5)))] = 100*min(x)
#x[sample(length(x), ceiling(0.5*length(x)^(1/5)))] = 100*max(x)
JHat1Med(x)
JHat1(x)
#c(JHat1(x), JHat1Med(x), JHat1MedBW(x), JHat1MedcT(x))
print(Sys.time()-start)

# Compare
start = Sys.time()
# all models are mean zero
J = 4 #AR(1) phi = 0.5 (sigma^2/(1-phi^2) * (1+2*(phi/(1-phi))))
#J = 2.25 #MA(1) theta = 0.5 (theta+1)^2 * sigma^2
#J = 9 #ARMA(1, 1) phi = 0.5 theta = 0.5  
nRep = 1000
n.all = c(100, 200, 400, 800)
out = array(NA, dim=c(nRep,length(n.all),2), 
            dimnames=list(paste0("iRep=",1:nRep),paste0("n=",n.all),c("EstMed","EstHat")))
for(iRep in 1:nRep){
  set.seed(iRep)
  #x.all = arima.sim(n=max(n.all), model=list(ar=.5))
  #x.all = arima.sim(n=max(n.all), model=list(ma=0.5))
  x.all = arima.sim(n=max(n.all), model=list(ar=0.5))

  for(i.n in 1:length(n.all)){
    n = n.all[i.n]
    x = x.all[1:n]
    x[sample(length(x), ceiling(0.5*length(x)^(1/3)))] = 3*min(x)
    x[sample(length(x), ceiling(0.5*length(x)^(1/3)))] = 3*max(x)
    #x[1:(length(x)^(2/3))] = 2*min(x)
    #ts.plot(x)
    out[iRep, i.n, 1] = JHat1Med(x)
    out[iRep, i.n, 2] = JHat1(x)
    
  }
  if(iRep%%100==0){
    cat(iRep," >> ")
  }
}
MSE = apply((out-J)^2,2:3,mean)
MSE
matplot(log(n.all), log(MSE), 
        type="l", pch=c("1","2"), col=c("red4","blue"),
        ylab="log(MSE)",xlab="log(n)",main="Comparison")
legend("topright", legend=c("Median","Average"), lty=1:2, col=c("red4","blue"))
print(Sys.time()-start)
