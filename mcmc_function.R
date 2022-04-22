

run_mcmc = function(model,       #model object to use
                    lab_data,    #in vitro dataset to use for model fitting
                    init.theta, #initial parameter values
                    proposal.sd = init.theta/100,  #standard deviation for parameter sampling
                    n.iterations = 1000,          #number of iterations
                    adapt.size.start = 200000,     #number of steps before adapting proposal distribution size
                    adapt.size.cooling = 0.99,     #cooling parameter for adaptation
                    adapt.shape.start = 200000,    #number of steps before adapting proposal distribution shape
                    verbose = FALSE){              
  
  target_function = function(theta){
    
    my_init.state = c(B = lab_data$B[1], A = lab_data$A[1])
    
    logval = dLogPosterior(fitmodel = model, theta = theta, init.state = my_init.state, 
                           data = lab_data, margLogLike = dTrajObs, log = TRUE)
    
    return(logval)
    
    
    
  }
  
  
  mcmc_fit = mcmcMH(target = target_function,
                    init.theta = init.theta,
                    proposal.sd = proposal.sd, 
                    n.iterations = n.iterations,
                    adapt.size.start = adapt.size.start,
                    adapt.size.cooling = adapt.size.cooling,
                    adapt.shape.start = adapt.shape.start,
                    verbose = verbose)
  
  mcmc_fit
  
}
