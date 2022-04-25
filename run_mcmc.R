
## SETUP ##########

if(!require(devtools, quietly = T)) install.packages("devtools")
if(!require(fitR, quietly = T)) devtools::install_github("sbfnk/fitR")
if(!require(ggplot2, quietly = T)) install.packages("ggplot2")
if(!require(scales, quietly = T)) install.packages("scales")
if(!require(here, quietly = T)) install.packages("here")
if(!require(dplyr, quietly = T)) install.packages("dplyr")
if(!require(reshape2, quietly = T)) install.packages("reshape2")

library(devtools)
library(fitR)
library(ggplot2)
library(scales)
library(here)
library(dplyr)
library(reshape2)

source(here::here("mcmc_function.R"))
model = readRDS(here::here("bac_model.rds"))


## RUN THE MODEL ##########

#parameter values
#mu: bacteria growth rate
#Nmax: carrying capacity
#beta: effect of antibiotic on bacteria
#gamma: antibiotic decay rate
theta = c(mu = 1.5, Nmax = 5e9, beta = 0.5, gamma = 0.3)

#times: run for 24h
times = seq(0, 24, 1)

#starting conditions
#10^4 cfu/mL bacteria, 10mg/L antibiotic
init.state = c(B = 1e4, A = 10)

#run model
results = as.data.frame(model$simulate(theta, init.state, times))

#plot number of bacteria over time
ggplot(results) +
  geom_line(aes(time, B)) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)),
                     limits = c(1, 5e10)) +
  theme_bw()


## GENERATE DUMMY LAB DATA ##########

#generate some dummy data by resampling model result
dummy_data = data.frame(time = results$time,
                        B1 = sapply(results$B, function(x) rpois(1, x)*sample(seq(0.5,1.5,0.1), 1)),
                        B2 = sapply(results$B, function(x) rpois(1, x)*sample(seq(0.5,1.5,0.1), 1)),
                        B3 = sapply(results$B, function(x) rpois(1, x)*sample(seq(0.5,1.5,0.1), 1)))

#plot dummy data
ggplot() +
  geom_line(data = results, aes(time, B)) +
  geom_point(data = dummy_data, aes(time, B1)) +
  geom_point(data = dummy_data, aes(time, B2)) +
  geom_point(data = dummy_data, aes(time, B3)) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)),
                     limits = c(1, 5e10)) +
  theme_bw()

#compute average bacterial concentration over time
#here we assume we know the antibiotic concentration is 10 mg/L, but we have no data
# to estimate the decay
lab_data = dummy_data
lab_data$B = rowMeans(lab_data[,-1])
lab_data = lab_data[,c(1,5)]
lab_data$A = 10

lab_data


## RUN MCMC ##########

#reminder, these are our "true" parameter values we're trying to find:
#theta = c(mu = 1.5, Nmax = 5e9, beta = 0.5, gamma = 0.3)

#pick some starting parameter values
init.theta = c(mu = 1, Nmax = 1e8, beta = 1, gamma = 0.8)

#run the model with these starting values
test_fit = model$simulate(theta = init.theta,
                          init.state = c(B = lab_data$B[1], A = lab_data$A[1]),
                          times = seq(0, 24, 1))

#compare model with starting values vs data
ggplot() +
  geom_line(data = lab_data, aes(time, B, colour = "Data")) +
  geom_line(data = test_fit, aes(time, B, colour = "Test fit")) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)),
                     limits = c(1, 5e10)) +
  theme_bw() +
  theme(legend.position = "bottom")

#choose number of iterations
n.iterations = 1000

#things you can try to play around with:
# 1) change the number of iterations
# 2) adjust proposal.sd (currently equal to 0.01 of each parameter initial value - what happens if you
#     divide by 10 instead? 1000?)
# 3) change the initial parameter values (init.theta)
# 4) try running adaptive mcmc (adjust the adapt.size.start to determine after how many iterations
#     the adaptive mcmc kicks in, adjust the adapt.shape start to determine after how many ACCEPTED
#     iterations the mcmc will consider it has reached the target distribution point and will focus
#     on exploring that target area instead of trying to move away)

#run the mcmc
mcmc_output = run_mcmc(model, lab_data,
                       init.theta = init.theta,
                       proposal.sd = init.theta/c(100,100,100,100),
                       n.iterations = n.iterations,
                       adapt.size.start = 2000000,
                       adapt.shape.start = 2000000)

#look at the trajectory of parameter values in the mcmc
mcmc_output$trace %>%
  as.data.frame() %>%
  select(-log.density) %>%
  melt() %>%
  ggplot() +
  facet_wrap(~variable, scales = "free") +
  geom_line(aes(x = c(1:(n.iterations*4)), y=value, group = variable)) +
  theme_bw() + 
  geom_hline(aes(yintercept = c(rep(theta[1], n.iterations),
                                rep(theta[2], n.iterations),
                                rep(theta[3], n.iterations),
                                rep(theta[4], n.iterations)),
                 color = "Target value")) + 
  geom_hline(aes(yintercept = c(rep(mean(mcmc_output$trace[,"mu"]), n.iterations),
                                rep(mean(mcmc_output$trace[,"Nmax"]), n.iterations),
                                rep(mean(mcmc_output$trace[,"beta"]), n.iterations),
                                rep(mean(mcmc_output$trace[,"gamma"]), n.iterations)),
                 color = "Average from mcmc")) +
  theme(legend.position = "bottom") +
  labs(x = "Step", y = "Value")

#look at the distribution of parameter values in the mcmc
mcmc_output$trace %>%
  as.data.frame() %>%
  select(-log.density) %>%
  melt() %>%
  ggplot() +
  geom_histogram(aes(value)) +
  facet_wrap(~variable, scales = "free") +
  theme_bw() +
  geom_vline(aes(xintercept = c(rep(theta[1], n.iterations),
                                rep(theta[2], n.iterations),
                                rep(theta[3], n.iterations),
                                rep(theta[4], n.iterations))),
             color = "blue") +
  labs(x = "Value", y = "Count")
  
#run the model with the best set of parameter values
best_fit = model$simulate(theta = mcmc_output$trace[which.max(mcmc_output$trace[,"log.density"]),],
                          init.state = c(B = lab_data$B[1], A = lab_data$A[1]),
                          times = seq(0, 24, 1))

#run the model with the average value of each parameter
best_fit_average = model$simulate(theta = apply(mcmc_output$trace, 2, mean),
                          init.state = c(B = lab_data$B[1], A = lab_data$A[1]),
                          times = seq(0, 24, 1))

#compare model with best fit vs model with average values vs dummy data
ggplot() +
  geom_line(data = lab_data, aes(time, B, colour = "Data")) +
  geom_line(data = best_fit, aes(time, B, colour = "Best fit")) +
  geom_line(data = best_fit_average, aes(time, B, colour = "Best fit (average)")) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)),
                     limits = c(1, 5e10)) +
  theme_bw() +
  theme(legend.position = "bottom")

