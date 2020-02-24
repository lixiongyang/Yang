# File name: rkc.R
# It uses to estimate the regression kink with an unknown constant threshold proposed by Hansen(2017)
# Hansen, B.E., 2017. Regression kink with an unknown threshold. J. Bus. Econ. Stat. 35 (2), 228–240.
# This code is based partly on the codes shared by Hansen(2017)

# Some useful functions: compute OLS estimates
reg <- function(X,y) {
  X <- qr(X)
  as.matrix(qr.coef(X,y))
}

pos.part <- function(x) x*(x>0)  #positive part
neg.part <- function(x) x*(x<0)  #negative part



rkc <- function(y,x,z, r0, r1, 
                   stp=0.1) {
  # model y = b1 (x-r)_- + b2 (x-r)_+ + b2*z + e
  # Input:
  #   y: outcome variable
  #   x: threshold dependent variable
  #   z: control variable
  #   stp:
  #   r0:
  #   r1:
  # gammas = seq(r0,r1,by=stp)	# Grid on Threshold parameter for estimation
  #   L.bt / U.bt: Lower and upper bounds for the threshold parameters gamma.
  #   L.dt / U.dt: Lower and upper bounds for dt
  #   L.gm / U.gm: Lower and upper bounds for gm
  #   tau1: Lower bound for the proportion of regime 1, i.e. (f'gm > 0)
  #   tau2: Upper bound for the proportion of regime 1, i.e. (f'gm > 0)
  #   eta: effective zero
  #   params: parameters for gurobi engine

  
  # Linear Model
  x0 = cbind(x,z)
  k0 = ncol(x0)
  kz = ncol(z)
  x00 = solve(crossprod(x0))  #inverse of x0'x0
  bols = x00%*%crossprod(x0,y)
  e0 = y - x0%*%bols
  sse0 = sum(e0^2)
  n = length(y)
  sigols = sse0/n
  v0 = x00%*%crossprod(x0*matrix(e0,n,k0))%*%x00*(n/(n-k0))
  seols = as.matrix(sqrt(diag(v0)))
  
  
  # Threshold Model
  gammas = seq(r0,r1,by=stp)	# Grid on Threshold parameter for estimation
  grid = length(gammas)
  
  sse = matrix(0,grid,1)
  k = kz + 3
  for (j in 1:grid) {
    gamj=gammas[j]
    x1 = cbind(neg.part(x-gamj),pos.part(x-gamj),z)
    e1 = y - x1%*%reg(x1,y)
    sse[j] = sum(e1^2)
  }
  gi = which.min(sse)
  gammahat = gammas[gi]
  ssemin = sse[gi]
  x1 = cbind(neg.part(x-gammahat),pos.part(x-gammahat),z)
  bt = reg(x1,y)
  et = y - x1%*% bt
  
  # compute se as in Theorem 2
  hg = - (x<gammahat)*bt[1] - (x>gammahat)*bt[2]
  x2 = cbind(x1,hg)
  hg2 = crossprod(cbind((x<gammahat),(x>gammahat)),et)
  xx2 = matrix(0,k,k)
  xx2[1:2,k]=hg2
  xx2[k,1:2]=t(hg2)
  xxi = solve(crossprod(x2) + xx2)
  v = xxi%*%crossprod(x2*matrix(et,n,k))%*%xxi*(n/(n-k))
  betahat = rbind(bt,gammahat)
  se = as.matrix(sqrt(diag(v)))
  sig = sum(et^2)/n
  wt = n*(sse0-ssemin)/ssemin
  wg = n*(sse-ssemin)/ssemin
  
  # print the estimates
  cat("----------------------------------------------","\n")
  cat("Linear Model, coefficients and error variance","\n")
  print(cbind(bols,seols),digits=2)
  cat("Error variance of linear model:","\n")
  print(sigols,digits=4)
  cat("----------------------------------------------","\n")
  
  
  cat("----------------------------------------------","\n")
  cat("Threshold Model Estimates, s.e.'s","\n")
  print(cbind(betahat,se),digits=3)
  cat("Error variance of kink model:","\n")
  print(sig,digits=4)
  cat("----------------------------------------------","\n")
  cat("F type statistic for kink effect:","\n")
  print(wt,digits = 4)
  
  return(list(bols=bols, seols=seols, sigols=sigols, betahat=bt, gammahat=gammahat, se=se, sig=sig, wt=wt))
}



##########################################################
# testing for kink effect
testkinkC <- function(y,x,z, boot=100, b0, bc, gammahat, w1, 
                      level=0.9) {
  
  # write a bootstrap for testing
  #b0 = rkc[1]$bols
  #bt = rkc[4]$betahat
  #gammahat = as.numeric(rkc[5])
  x1 = cbind(neg.part(x-gammahat),pos.part(x-gammahat),z)
  ec = y - x1%*% bc
  n = length(y)  #lag
 
  
  
  Fb = matrix(0,boot)
  xli = cbind(x,z)
  
  for (j in 1:boot) {
    yb = xli%*%b0 + ec*rnorm(n)
    source("rkc.R")
    rkcb <- rkc(yb,x,z,r0=r0,r1=r1,stp=stp)
    sse0b = as.numeric(rkcb[3])
    ssekb = as.numeric(rkcb[7])
    Fb[j]=n*(sse0b-ssekb)/ssekb
  }
  
  pv = mean(Fb > matrix(w1,boot,1))
  crit = quantile(Fb,probs=0.9)
  
  cat('critical value   =', crit, '\n')
  cat("----------------------------------------------","\n")
  cat("statistic(linear vs constant kink)=", "\n")
  print(w1)
  cat('p-value   =', pv, '\n')
  cat("----------------------------------------------","\n")
  return(list(pv, crit))
}