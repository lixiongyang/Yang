# Load in Data
rm(list=ls())
setwd("E:\\R-econometrics\\cthresh-Hansen-2017")
growth <- read.table("usdata.txt",header=TRUE)
n = nrow(growth)
year = growth[2:n,1]
gdp  = growth[2:n,3]
gdp1 = growth[1:(n-1),3]
debt1= growth[1:(n-1),2]

# Time-Series Data Plots, Figure 1ab in paper
#windows()
plot(year,gdp,type="l",ylab="GDP Growth Rate")
#savePlot(file="fig1a.eps",type="eps",dev.cur())

#windows()
plot(year,debt1,type="l",ylab="Debt/GDP")
#savePlot(file="fig1b.eps",type="eps",dev.cur())

# Define variables
y = gdp
x = debt1
n = length(y)
z = cbind(gdp1,matrix(1,n,1))    # controls and constant


#load functions
setwd("E:\\R-econometrics\\MyRegressionKink")
source("rkc.R")
# set grid search paramaters
r0 = 10
r1 = 71
stp = 0.1
# estimate the model
rkc <- rkc(y,x,z,r0=r0,r1=r1,stp=stp)  
#wt2 = n*(sse0-ssemin)/ssemin
sse0 = as.numeric(rkc[3])
ssek = as.numeric(rkc[7])
wt = 218*(sse0-ssek)/ssek


# set parameters a bootstrap for testing kink effect
b0 = rkc[1]$bols   # coefficients of linear model: null model
bc = rkc[4]$betahat  # coefficients of kink model proposed by Hansen(2017)
gammahat = as.numeric(rkc[5])
w1 = rkc[8]   # statistic for kink effect W1
level = 0.9 # significant level
boot = 10

testkinkC <- testkinkC(y,x,z, boot=boot, b0, bc, gammahat, w1, level=level)
