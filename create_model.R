
model_name = "bac model"
model_state.names = c("B","A")
model_theta.names = c("mu", "Nmax", "beta", "gamma")

model_simulateDeterministic <- function(theta, init.state, times) {
  
  model_ode = function(times, y, parms){
    
    with(as.list(c(y, parms)), {
      
      dBdt = mu*(1 - B/Nmax)*B - beta*A*B
      
      dAdt = -gamma*A
      
      return(list(c(dBdt, dAdt)))
      
    })
  }
  
  trajectory = data.frame(ode(y = init.state,
                              times = times,
                              func = model_ode,
                              parms = theta,
                              method = "ode45"))
  
  return(trajectory)
  
}


## function to compute log-prior
model_prior <- function(theta, log = TRUE) {
  
  log.prior.mu <- dunif(theta[["mu"]], min = 0.01, max = 10, log = TRUE)
  log.prior.Nmax <- dunif(theta[["Nmax"]], min = 1e7, max = 1e11, log = TRUE)
  log.prior.beta <- dunif(theta[["beta"]], min = 0.01, max = 10, log = TRUE)
  log.prior.gamma <- dunif(theta[["gamma"]], min = 0.01, max = 10, log = TRUE)
  
  log.sum <- log.prior.mu + log.prior.Nmax + log.prior.beta + log.prior.gamma
  
  return(ifelse(log, log.sum, exp(log.sum)))
}



## function to compute the likelihood of one data point
model_pointLike <- function(data.point, model.point, theta, log = FALSE){
  
  dpoisB = dpois(x = round(data.point[["B"]]/(10^(max(nchar(as.character(round(model.point[["B"]]))),2)-2))),
                 lambda = model.point[["B"]]/(10^(max(nchar(as.character(round(model.point[["B"]]))),2)-2)),
                 log = log)
  
  ## the prevalence is observed through a Poisson process
  return(dpoisB)
}

## create deterministic SIR fitmodel
model <- fitR::fitmodel(
  name = model_name,
  state.names = model_state.names,
  theta.names = model_theta.names,
  simulate = model_simulateDeterministic,
  dprior = model_prior,
  dPointObs = model_pointLike)

saveRDS(model, here::here("bac_model.rds"))

