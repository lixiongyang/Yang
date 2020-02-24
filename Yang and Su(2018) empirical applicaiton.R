
# Load in Data
rm(list=ls())
setwd("E:\\R-econometrics\\cthresh-Hansen-2017")
growth <- read.table("usdata.txt",header=TRUE)
n = nrow(growth)
year = growth[2:n,1]
gdp  = growth[2:n,3]
gdp1 = growth[1:(n-1),3]
debt1= growth[1:(n-1),2]
q = growth[2:n,4]
#q = gdp1

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
source("rkt.R")
# set grid search paramaters
r01 = 10
r02 = 71
stp1 = 1
r11 = -5
r12 = 5
stp2 = 0.1
# estimate the model
rkt <- rkt(y,x,z,q,r01=r01,r02=r02,r11=r11,r12=r12,stp1=stp1,stp2=stp2)  


#wt2 = n*(sse0-ssemin)/ssemin
sse0 = as.numeric(rkt[3])
ssek = as.numeric(rkt[7])
wt = n*(sse0-ssek)/ssek


# write a bootstrap for testing kink effect
b0 = rkt[1]$bols
bt = rkt[4]$bt   #coefficients for kink model with a varying threshold
gammahat0 = as.numeric(rkt[5])
gammahat1 = as.numeric(rkt[6])
w2 = rkt[8]  # statistic for kink effect (state-dependent threshold)
boot = 10 # boostrap replications

# testing for kink effect based on Yang and Su (2018, JIMF)
testkinkT <- testkinkT(y,x,z,q, boot=boot, b0, bt, gammahat0, gammahat1, 
                      w2, level=0.9)


