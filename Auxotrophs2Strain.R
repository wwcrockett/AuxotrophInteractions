# Load the ODE Solver Package!
require(deSolve)

# 2 Strain ODE Model:
ode.Model = function(time, states, params){
  
  # Extract the state variables
  x1 = states[1]     # Strain 1 Population
  x2 = states[2]     # Strain 2 Population
  
  # Extract the parameter values
  beta = params[1]    # buffer beta for low populations
  k = params[2]    # total population carrying capacity k
  c12 = params[3]    # effect of x2 on x1
  c21 = params[4]    # effect of x1 on x2
  
  # Evaluate the differential equation at time t

  dx1dt = c12*x2*(x1/(x1+beta))*(1-(x1+x2)/k)   # Equation for Strain 1
  dx2dt = c21*x1*(x2/(x2+beta))*(1-(x1+x2)/k)   # Equation for Strain 2

  # Combine into a single vector
  dXdt = c(dx1dt,dx2dt)
  
  # Return the calculations as a list as required by ODE solver
  return(list(dXdt))
  
}
# Example of the model ----

# Set up time vector
time = seq(0,5500,by=1)

# Set up vector of parameter values
beta = 1
k = 1e9
c12 = 0.005
c21 = 0.004

params = c(beta,k,c12,c21)      # Parameter Vector

# Set initial values for state variables
x1o = 1e7
x2o = 1e7

states0 = c(x1o,x2o)      # State Variable Vector

# Numerically solve the ODE:

numSolutions = ode(states0,time,ode.Model,params)

# Plot output from model:
matplot(numSolutions[1:3000,1], numSolutions[1:3000,2:3], type=c("l","l"), 
        lty=c(1,1), lwd=c(2,2), xlab="",ylab="",col=c("blue","orange"), log = 'y')

# Add a legend:
legend('right',legend=c("Strain 1","Strain 2"), col=c("blue","orange"), 
       lty=c(1,1), lwd=c(2,2), cex=0.75)

# Add a title and axes labels
title("Population Dynamics",xlab="time",ylab="Population Density [Cells/mL]")



# Function for evaluating the sum of squared errors (SSE) ----
# between model output and data
SSE = function(parmstofit,expdata){
  
  # Exponentiate parameters to undo log (and make them positive)
  parmstofit = exp(parmstofit)
  # We are not fitting initial population, beta, or K:
  x1o = 1e7
  x2o = 1e7
  states0 = c(x1o,x2o)      # State Variable Vector
  
  beta = 1
  k = 1e9
  
  # We are fitting the interaction coefficients:
  c12 = parmstofit[1]
  c21 = parmstofit[2]
  
  params = c(beta,k,c12,c21)      # Parameter Vector
  # Specify output times to coincide with those of data
  time = seq(0,6000,by=1)
  # Find solution to ODE
  out = ode(states0,time,ode.Model,params)
  # Calculate sum of squared errors between 
  # Fold growth for each population
  fold1 = out[5500,2]/out[1,2]
  fold2 = out[5500,3]/out[1,3]
  sumsqerrors = (fold1 - expdata[1])^2+(fold2-expdata[2])^2
  return(sumsqerrors)
}

# Load data from text file into dataframe ----
fold.growths = read.csv("foldgrowth_2member.csv")
fold.growths = fold.growths[1:91,1:7]


fits1 = matrix(0,nrow = dim(fold.growths)[1], ncol = 6)
fits2 = matrix(0,nrow = dim(fold.growths)[1], ncol = 6)
Pop = matrix(0,nrow = dim(fold.growths)[1], ncol= 1)
colnames(fits1) <- c("c12fit","c12paper","c12dif","fgfit","fgpaper","fgdif")
rownames(fits1) <- fold.growths[,2]
colnames(fits2) <- c("c21fit","c21paper","c21dif","fgfit","fgpaper","fgdif")
rownames(fits2) <- fold.growths[,5]
colnames(Pop) <- "FG"
rownames(Pop) <- fold.growths[,1]

beta = 1
K = 1e9
time = seq(0,5500,by=1)

# Fit to Fold Growth for each Combination of strains: ----
for (i in fold.growths[,1]){
  
  # Initial guesses for parameters
  c12 = 1e-7
  c21 = 1e-7

  params = c(c12,c21)      # Parameter Vector

  # Log transform so that parameters are unconstrained 
  # (can be negative) when passed to optim
  init.guess = log(params)
  #init.guess = params
  # Using R's optimization algorithm to minimize SSE 
  modelfit = optim(init.guess,SSE,expdata=c(fold.growths[i,3],fold.growths[i,6]),method = "Nelder-Mead") 
  c12 = exp(modelfit$par)[1]
  c21 = exp(modelfit$par)[2]
  fits1[i,"c12paper"] = fold.growths[i,4]
  fits2[i,"c21paper"] = fold.growths[i,7]
  fits1[i,"c12fit"] = c12
  fits2[i,"c21fit"] = c21
  fits1[i,"c12dif"] = (c12-fits1[i,"c12paper"])/fits1[i,"c12paper"]
  fits2[i,"c21dif"] = (c21-fits2[i,"c21paper"])/fits2[i,"c21paper"]
  ps = c(beta,K,c12,c21)
  numSolutions = ode(states0,time,ode.Model,ps)
  fits1[i,"fgfit"] = numSolutions[5500,2]/numSolutions[1,2]
  fits2[i,"fgfit"] = numSolutions[5500,3]/numSolutions[1,3]
  Pop[i,"FG"] = (numSolutions[5500,3] + numSolutions[5500,2])/(numSolutions[1,2]+numSolutions[1,3])
  fits1[i,"fgpaper"] = fold.growths[i,3]
  fits2[i,"fgpaper"] = fold.growths[i,6]
  fits1[i,"fgdif"] = fits1[i,"fgfit"] - fits1[i,"fgpaper"]
  fits2[i,"fgdif"] = fits2[i,"fgfit"] - fits2[i,"fgpaper"]
}

# Save Fits ----
df = data.frame(row.names(fits1),row.names(fits2),fits1,fits2)
rownames(df) <- NULL
write.csv(df,"Bestfits.csv")

# Plot comparisons to paper's fits ----

plot(fits1[,"c12paper"],fits1[,"c12dif"],xlab="c12 Paper Fits",ylab = "Difference in Fit")


plot(fits1[,"fgpaper"],fits1[,"fgdif"], xlab = "Measured Strain 1 Fold Growth",ylab="Strain 1 Fold Growth Difference")
title("Difference in Strain 1 Fold Growths")


plot(fits2[,"c21paper"],fits2[,"c21dif"],xlab="c21 Paper Fits",ylab = "Difference in Fit")


plot(fits2[,"fgpaper"],fits2[,"fgdif"], xlab = "Measured Strain 2 Fold Growth",ylab="Strain 2 Fold Growth Difference")
title("Difference in Strain 2 Fold Growths")


