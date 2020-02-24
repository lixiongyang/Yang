
rm(list=ls())
#load functions
setwd("E:\\R-econometrics\\MyRegressionKink")
source("rkt.R")

# generate data

n=200
x = rnorm(n)
q = rnorm(n)
rt = 0.2 - 0.5*q
z = rnorm(n)
x1 = cbind(neg.part(x-rt),pos.part(x-rt),z)
b0 =c(1,2,1)
y = x1%*%b0 + rnorm(n)



# set grid search paramaters
r01 = 0
r02 = 2
stp1 = 0.1
r11 = -10
r12 = 5
stp2 = 0.1
# estimate the model with a state-dependent threshold
rkt <- rkt(y,x,z,q,r01=r01,r02=r02,r11=r11,r12=r12,stp1=stp1,stp2=stp2)  