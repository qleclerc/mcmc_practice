
#how you put mcmc in practice with a model/story
#setup small model
#produce random data
#try to refit model to data
#small run first
#then long run

#further considerations to show complexities must think about (parameters correlated, data processing,
#fitting on multiple curves at once...)

library(deSolve)
library(ggplot2)
library(scales)

bac_model = function(times, y, parms){

  with(as.list(c(y, parms)), {
  
  dBdt = mu*(1 - B/Nmax)*B - beta*A*B
  
  dAdt = -gamma*A
  
  return(list(c(dBdt, dAdt)))
  
  })
}

parms = c(mu = 1.5, Nmax = 1e9, beta = 0.5, gamma = 0.3)

times = seq(0, 24, 1)

yinit = c(B = 1e4, A = 10)

results = as.data.frame(ode(func = bac_model, times = times, y = yinit, parms = parms))

ggplot(results) +
  geom_line(aes(time, B)) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)),
                     limits = c(1, 5e9)) +
  theme_bw()


dummy_data = data.frame(time = results$time,
                        B1 = sapply(results$B, function(x) rpois(1, x)*sample(seq(0.6,1.4,0.1), 1)),
                        B2 = sapply(results$B, function(x) rpois(1, x)*sample(seq(0.6,1.4,0.1), 1)),
                        B3 = sapply(results$B, function(x) rpois(1, x)*sample(seq(0.6,1.4,0.1), 1)))

ggplot() +
  geom_line(data = results, aes(time, B)) +
  geom_point(data = dummy_data, aes(time, B1)) +
  geom_point(data = dummy_data, aes(time, B2)) +
  geom_point(data = dummy_data, aes(time, B3)) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)),
                     limits = c(1, 5e9)) +
  theme_bw()


