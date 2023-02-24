
library(BayesianTools)
library(dplyr)
library(deSolve)
library(ggplot2)

#this wrapper should run your model for a given set of parameters, starting conditions, and times,
# and return a trajectory
#for greater efficiency, you likely want to define the function in environment, then just
# call it in the wrapper, instead of defining it in the wrapper
SIR_wrapper = function(parameters, init.state, times){      
  
  #replace this with your model
  SIR_model = function(t, pop, param) {
    S = pop[1]
    I = pop[2]
    R = pop[3]
    
    B=param[1]
    G=param[2]
    u0=param[3]
    u1=param[4]
    
    N=S+I+R
    
    dS = -B*S*I/N + u0*N - u1*S
    dI = B*S*I/N  - G*I - u1*I
    dR =  G*I - u1*R
    
    res=c(dS, dI, dR)
    
    list(res)
  }
  
  trajectory = data.frame(lsoda(y = init.state,
                                times = times,
                                func = SIR_model,
                                parms = parameters))
  
  return(trajectory)
}

#generate some observation data to test the fitting
beta=3
gamma=0.5
u0=0.01
u1=0.01

I0=1
S0=1000-10
R0=0

dt=1
Tmax=30

times = seq(from=0,to=Tmax,by=dt)
init.state = c(S=S0,I=I0,R=R0) 
parameters = c(beta,gamma,u0,u1)

obs = SIR_wrapper(parameters = parameters, init.state = init.state, times = times)

ggplot(obs) +
  geom_line(aes(times, S, colour = "S")) +
  geom_line(aes(times, I, colour = "I")) +
  geom_line(aes(times, R, colour = "R")) +
  theme_bw() +
  labs(x = "Time (days)", y = "Individuals", colour = "")

#add some noise to the observation data
obs[,-1] = apply(obs[,-1], c(1,2), function(x) rpois(1,x))

ggplot(obs) +
  geom_line(aes(times, S, colour = "S")) +
  geom_line(aes(times, I, colour = "I")) +
  geom_line(aes(times, R, colour = "R")) +
  theme_bw() +
  labs(x = "Time (days)", y = "Individuals", colour = "")


#this is basically just a dataframe with all your parameters, alongside the 
# lower and upper values to use in the MCMC below
refPars = data.frame(best = c(5, 0.5, 0.01, 0.01),
                     lower = c(1, 0.5, 0.01, 0.01),
                     upper = c(10, 1, 0.5, 0.01))
rownames(refPars) = c("beta", "gamma", "u0", "u1")

#the advantage of this structure is that you can choose which parameter to estimate
#here, we only want to estimate beta (parameter 1). The sampler will then just pick
# the "best" values we saved in refPars for all other parameters
parSel = c(1)

#here is the likelihood function
#this should take as input the parameters which the sampler will propose, run the 
# model, calculate a likelihood value, and return this point estimate of likelihood
likelihood = function(par){
  
  #parameters that we are not estimating are set on default values defined in refPars
  x = refPars$best
  x[parSel] = par
  
  #run the model wrapper with the parameters to test
  predicted = SIR_wrapper(parameters = x,
                          init.state = c(S = obs$S[1], I = obs$I[1], R = obs$R[1]),
                          times = obs$time)
  
  #quick fix to avoid a likelihood crash with Poisson if some values are < 0
  predicted[predicted<0] = 0
  
  #calculate likelihood
  llValues = dpois(x = obs$I,
                   lambda = predicted$I,
                   log = T)
  
  #return a point estimate of likelihood
  return(sum(llValues))
}

#creating the prior(s) for the parameter(s) you want to estimate, based on what you stored
# in refPars
prior = createUniformPrior(lower = refPars$lower[parSel], 
                           upper = refPars$upper[parSel],
                           best = refPars$best[parSel])

bayesianSetup = createBayesianSetup(likelihood, prior, names = rownames(refPars)[parSel])

#settings for the sampler, see documentation for more info
settings = list(iterations = 10000)

#run the mcmc
out = runMCMC(bayesianSetup = bayesianSetup, settings = settings)

#plot posterior
plot(out)

#summary of mcmc run
summary(out)

#a quick check to compare the initial guess output, predicted output, and observed output
#initial guess output:
beta = refPars["beta","best"] #out initial guess was stored in refPars
parameters = c(beta,gamma,u0,u1) #redefine parameter vector
#run model with estimated beta
init_guess = SIR_wrapper(parameters = parameters, init.state = init.state, times = times)


#predicted output:
beta = median(getSample(out)) #take posterior median for beta
parameters = c(beta,gamma,u0,u1) #redefine parameter vector
#run model with estimated beta
pred = SIR_wrapper(parameters = parameters, init.state = init.state, times = times)

#compare predicted I and observed I
ggplot() +
  geom_line(data = init_guess, aes(times, I, colour = "initial_guess_I")) +
  geom_line(data = pred, aes(times, I, colour = "predicted_I")) +
  geom_line(data = obs, aes(times, I, colour = "observed_I")) +
  theme_bw() +
  labs(x = "Time (days)", y = "Individuals", colour = "")
