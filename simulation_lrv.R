# generate data from AR(1)
x = arima.sim(n=1000, model=list(ar=.5))
ts.plot(x)

# autocovariance
cHat = c(acf(x, type="cov")$acf)

# long-run variance 
JHat = function(x, b=NA){
	if(is.na(b)){
		b = ceiling(2*length(x)^(1/3))
	}
	cHat = c(acf(x, type="cov", lag.max=b, plot=FALSE)$acf)	
	kernel = 1-(0:b)/b
	temp = cHat*kernel
	temp[1] + 2*sum(temp[-1])
}
x = arima.sim(n=1000, model=list(ar=.5))
J = 4
JHat(x)

kBT = function(x, T){
  if(abs(x) <= 1)
    return(1-abs(x))
  else
    return(0)
}
JCla2 = function(x, b=NA){
  T = length(x)
  if(is.na(b)){
    b = ceiling(2*length(x)^(1/3))
  }
  
  res = 0
  for(k in 0:(T-1)){
    gamma = sum(x[(k+1):T]*x[1:(T-k)])
    #if(k <= 20)
      #cat(k, ": acf", kBT(k/b, T)*gamma/T, "\n")
    if(k == 0)
      res = res + kBT(k/b, T)*gamma/T
    else
      res = res + 2*kBT(k/b, T)*gamma/T
    #if(k %% 100 == 0)
    #cat(k, ": res =", res, "\n")
  }
  res
}

# replicate for 1000 times 
nRep = 1000
n.all = c(100,1000,10000)
out = array(NA, dim=c(nRep,length(n.all),2), 
			dimnames=list(paste0("iRep=",1:nRep),paste0("n=",n.all),c("Est1","Est2")))
for(iRep in 1:nRep){
	set.seed(iRep)
	x.all = arima.sim(n=max(n.all), model=list(ar=.5))
	for(i.n in 1:length(n.all)){
		n = n.all[i.n]
		x = x.all[1:n]
		out[iRep, i.n, 1] = JHat(x)
		out[iRep, i.n, 2] = JCla2(x)
	}
	if(iRep%%100==0){
		cat(iRep," >> ")
	}
}
out
MSE = apply((out-J)^2,2:3,mean)
MSE
matplot(log(n.all), log(MSE), 
	type="b", pch=c("1","2"), col=c("red4","green3"),
	ylab="log(MSE)",xlab="log(n)",main="Comparison")




