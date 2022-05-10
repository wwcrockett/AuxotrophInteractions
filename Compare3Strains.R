#### Auxotrophs 3 Strain:

# Using the Interactions found by fitting the 2 strain population to model 
# 3 strain populations.

# Load the ODE Solver Package!
require(deSolve)
# Read in interaction coefficients produced by Auxotrophs2Strain.R:----
Fits = read.csv("BestFits.csv")

# Defining Limiting Model and Lotka Volterra Model ----
# Also tried Multiplying model, but wasn't very useful
ode.Limiting = function(time, states, params){
  
  # Extract the state variables
  x1 = states[1]     # Strain 1 Population
  x2 = states[2]     # Strain 2 Population
  x3 = states[3]
  # Extract the parameter values
  beta= params[1]    # buffer beta for low populations
  k   = params[2]    # total population carrying capacity k
  c12 = params[3]    # effect of x2 on x1
  c21 = params[4]    # effect of x1 on x2
  c13 = params[5]
  c31 = params[6]
  c23 = params[7]
  c32 = params[8]
  # Evaluate the differential equation at time t
  dx1dt = min(c12*x2,c13*x3)*(x1/(x1+beta))*(1-(x1+x2+x3)/k)   # Equation for Strain 1
  dx2dt = min(c21*x1,c23*x3)*(x2/(x2+beta))*(1-(x1+x2+x3)/k)   # Equation for Strain 2
  dx3dt = min(c31*x1,c32*x2)*(x3/(x3+beta))*(1-(x1+x2+x3)/k)   # Equation for Strain 3
  # Combine into a single vector
  dXdt = c(dx1dt,dx2dt,dx3dt)
  
  # Return the calculations as a list as required by ODE solver
  return(list(dXdt))
}

ode.GenLV = function(time, states, params){
  
  # Extract the state variables
  x1 = states[1]     # Strain 1 Population
  x2 = states[2]     # Strain 2 Population
  x3 = states[3]
  # Extract the parameter values
  beta= params[1]    # buffer beta for low populations
  k   = params[2]    # total population carrying capacity k
  c12 = params[3]    # effect of x2 on x1
  c21 = params[4]    # effect of x1 on x2
  c13 = params[5]
  c31 = params[6]
  c23 = params[7]
  c32 = params[8]
  # Evaluate the differential equation at time t
  dx1dt = (c12*x2+c13*x3)*(x1/(x1+beta))*(1-(x1+x2+x3)/k)   # Equation for Strain 1
  dx2dt = (c21*x1+c23*x3)*(x2/(x2+beta))*(1-(x1+x2+x3)/k)   # Equation for Strain 2
  dx3dt = (c31*x1+c32*x2)*(x3/(x3+beta))*(1-(x1+x2+x3)/k)   # Equation for Strain 3
  # Combine into a single vector
  dXdt = c(dx1dt,dx2dt,dx3dt)
  
  # Return the calculations as a list as required by ODE solver
  return(list(dXdt))
}

ode.Mult = function(time, states, params){
  
  # Extract the state variables
  x1 = states[1]     # Strain 1 Population
  x2 = states[2]     # Strain 2 Population
  x3 = states[3]
  # Extract the parameter values
  beta= params[1]    # buffer beta for low populations
  k   = params[2]    # total population carrying capacity k
  c12 = params[3]    # effect of x2 on x1
  c21 = params[4]    # effect of x1 on x2
  c13 = params[5]
  c31 = params[6]
  c23 = params[7]
  c32 = params[8]
  # Evaluate the differential equation at time t
  dx1dt = (c12*x2*c13*x3)*(x1/(x1+beta))*(1-(x1+x2+x3)/k)   # Equation for Strain 1
  dx2dt = (c21*x1*c23*x3)*(x2/(x2+beta))*(1-(x1+x2+x3)/k)   # Equation for Strain 2
  dx3dt = (c31*x1*c32*x2)*(x3/(x3+beta))*(1-(x1+x2+x3)/k)   # Equation for Strain 3
  # Combine into a single vector
  dXdt = c(dx1dt,dx2dt,dx3dt)
  
  # Return the calculations as a list as required by ODE solver
  return(list(dXdt))
}

# Generate all unique combinations of 3 strains:  -----------------------------
AAs = unique(sort(append(Fits[,"row.names.fits1."],Fits[,"row.names.fits2."])))

Triplets = matrix(0,nrow = 364,ncol = 3)
number =0
for (i in 1:(length(AAs)-2)){
  for(j in (i+1):(length(AAs)-1)){
    for(k in (j+1):length(AAs)){
      number = number+1
      Triplets[number,1] = AAs[i]
      Triplets[number,2] = AAs[j]
      Triplets[number,3] = AAs[k]
    }
  }
}

# Read in Measured Fold Growth for 3 strain populations:
Measured = read.csv("foldgrowth_3member.csv")

# Finding the fits and simulating ----

foldtotLV = rep(0,number)
foldtotLim = rep(0,number)
foldtotMult = rep(0,number)
foldtotMes = rep(0,number)
LimDif = rep(0,number)
LVDif = rep(0,number)
MultDif = rep(0,number)
# Set up time vector
time = seq(0,5500,by=1)

# Set up vector of parameter values
beta = 1
k = 1e9

# Set initial values for state variables
x1o = 1e7
x2o = 1e7
x3o = 1e7

states0 = c(x1o,x2o,x3o)      # State Variable Vector

for (i in 1:number){
  
  # Find the pair interaction coefficients:
  rows12 = Fits[Fits$row.names.fits1.==Triplets[i,1],]
  row12 = rows12[rows12$row.names.fits2.==Triplets[i,2],]
  c12 = as.numeric(row12["c12fit"])
  c21 = as.numeric(row12["c21fit"])
  if (length(row12$row.names.fits1.)==0){
    rows12 = Fits[Fits$row.names.fits1.==Triplets[i,2],]
    row12 = rows12[rows12$row.names.fits2.==Triplets[i,1],]
    c12 = as.numeric(row12["c21fit"])
    c21 = as.numeric(row12["c12fit"])
  }
  rows13 = Fits[Fits$row.names.fits1.==Triplets[i,1],]
  row13 = rows12[rows13$row.names.fits2.==Triplets[i,3],]
  c13 = as.numeric(row12["c12fit"])
  c31 = as.numeric(row12["c21fit"])
  if (length(row13$row.names.fits1.)==0){
    rows13 = Fits[Fits$row.names.fits1.==Triplets[i,3],]
    row13 = rows12[rows13$row.names.fits2.==Triplets[i,1],]
    c13 = as.numeric(row12["c21fit"])
    c31 = as.numeric(row12["c12fit"])
  }
  rows23 = Fits[Fits$row.names.fits1.==Triplets[i,2],]
  row23 = rows12[rows13$row.names.fits2.==Triplets[i,3],]
  c23 = as.numeric(row12["c12fit"])
  c32 = as.numeric(row12["c21fit"])
  if (length(row23$row.names.fits1.)==0){
    rows23 = Fits[Fits$row.names.fits1.==Triplets[i,3],]
    row23 = rows12[rows13$row.names.fits2.==Triplets[i,2],]
    c23 = as.numeric(row12["c21fit"])
    c32 = as.numeric(row12["c12fit"])
  }
  params = c(beta,k,c12,c21,c13,c31,c23,c32)
  
  # Find the Measured Total Fold Growth:
  # Dataset is not very organized, have to look through each strain combination ----
  rows1 = Measured[Measured$Strain.1==Triplets[i,1],]
  rows2 = rows1[rows1$Strain.2==Triplets[i,2],]
  rows3 = rows2[rows2$Strain.3==Triplets[i,3],]
  if (length(rows3$Strain.1)== 0){
    rows2 = rows1[rows1$Strain.2==Triplets[i,3],]
    rows3 = rows2[rows2$Strain.3==Triplets[i,2],]
  }
  if (length(rows3$Strain.1)== 0) {
    rows1 = Measured[Measured$Strain.1==Triplets[i,2],]
    rows2 = rows1[rows1$Strain.2==Triplets[i,1],]
    rows3 = rows2[rows2$Strain.3==Triplets[i,3],]
  }
  if (length(rows3$Strain.1)== 0) {
    rows1 = Measured[Measured$Strain.1==Triplets[i,2],]
    rows2 = rows1[rows1$Strain.2==Triplets[i,3],]
    rows3 = rows2[rows2$Strain.3==Triplets[i,1],]
  }
  if (length(rows3$Strain.1)== 0) {
    rows1 = Measured[Measured$Strain.1==Triplets[i,3],]
    rows2 = rows1[rows1$Strain.2==Triplets[i,1],]
    rows3 = rows2[rows2$Strain.3==Triplets[i,2],]
  }
  if (length(rows3$Strain.1)== 0) {
    rows1 = Measured[Measured$Strain.1==Triplets[i,3],]
    rows2 = rows1[rows1$Strain.2==Triplets[i,2],]
    rows3 = rows2[rows2$Strain.3==Triplets[i,1],]
  }
  foldtotMes[i] = rows3$Strain.1.2.3.fold.growth..T84. # Total measured fold growth
  # Use pair coefficients for limiting model: ----
  LimSolutions = ode(states0,time,ode.Limiting,params)
  foldtotLim[i] = sum(LimSolutions[5500,2:4])/sum(LimSolutions[1,2:4])
  LimDif[i] = foldtotMes[i] - foldtotLim[i]
  
  # Use pair coefficients for generalized Lotka-Volterra model:
  LVSolutions = ode(states0,time,ode.GenLV,params)
  foldtotLV[i] = sum(LVSolutions[5500,2:4])/sum(LVSolutions[1,2:4])
  LVDif[i] = foldtotMes[i] - foldtotLV[i]
  
  # Use pair coefficients for multiplicative model:
  MultSolutions = ode(states0,time,ode.Mult,params)
  foldtotMult[i] = sum(MultSolutions[5500,2:4])/sum(MultSolutions[1,2:4])
  MultDif[i] = foldtotMes[i] - foldtotMult[i]
}
# Plot Total Fold-Growths for each Model----
plot(foldtotMes,LimDif, xlab = "Total Measured Fold-Growth",ylab="Fold-Growth Difference")
title("Limiting Model Difference in Fold-Growths")
plot(foldtotMes,LVDif, xlab = "Total Measured Fold-Growth",ylab="Fold-Growth Difference")
title("Lotka-Volterra Model Difference in Fold-Growths")
plot(foldtotMes,MultDif)
scatter.smooth(foldtotMes,LimDif)

LimDifMean = sum(LimDif)/number
LimMeans = rep(LimDifMean,number)
LimDifSigma = sqrt(sum((LimDif-LimDifMean)^2)/(number-1))
LVDifMean = sum(LVDif)/number
LVDifSigma = sqrt(sum((LVDif-LVDifMean)^2)/(number-1))
