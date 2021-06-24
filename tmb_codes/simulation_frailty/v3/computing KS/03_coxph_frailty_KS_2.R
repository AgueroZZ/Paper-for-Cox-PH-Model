# ## Coxph regression with frailty ###
# lib_loc <- '/home/ziang/lib'
# library(tidyverse, lib = lib_loc)
# library(aghq, lib = lib_loc)
# library(mgcv)
# library(Matrix)
# library(rstan, lib = lib_loc)
# library(TMB, lib = lib_loc)
# library(INLA, lib = lib_loc)
# library(tmbstan, lib = lib_loc)
# library(foreach, lib = lib_loc)
# library(doMC, lib = lib_loc)
# library(parallel)
# library(foreach)
# library(abcoxph, lib = lib_loc)
# library(mvQuad, lib = lib_loc)
# 
library(tidyverse)
library(aghq)
library(mgcv)
library(Matrix)
library(rstan)
library(TMB)
library(INLA)
library(tmbstan)
library(foreach)
library(doMC)
library(parallel)
library(foreach)
library(abcoxph)
library(mvQuad)

precompile()
TEXT_SIZE = 25
ncores = 10
registerDoMC(ncores)


### Simulating function:
K = 60
beta = 0.2
sd = 0.8
N <- 2

Simulate_baseline3 <- function(timelim = 300, breaks = 0.001, cut = 5){
  timelim <- timelim
  tdom <- seq(0, timelim, by = breaks)
  haz <- rep(0, length(tdom))
  cut <- cut
  for (i in 1:cut) {
    low <- as.numeric(quantile(tdom,(i-1)/cut))
    high <- as.numeric(quantile(tdom,(i)/cut))
    if(i %% 2 == 1){
      a <- runif(1,0,1)
      if(a > 0.3) haz[tdom<=high & tdom > low] <- 0.1
      else {
        c <- tdom[tdom<=high & tdom > low]
        haz[tdom<=high & tdom > low] <-0.01
      }
    }
    if(i %% 2 == 0){
      a <- runif(1,0,1)
      if(a > 0.8){
        c <- tdom[tdom<=high & tdom > low]
        haz[tdom<=high & tdom > low] <- 0.25
      }
      else{
        haz[tdom<=high & tdom > low] <- sample(c(0.05,0.15),size = 1,prob = c(0.5,0.5))
      }
    }
  }
  baseline <- data.frame(time = tdom, hazard = haz, timelim = timelim)
}
### Try this baseline on the smoothing example as well
set.seed(1234)
baseline <- Simulate_baseline3()

Simulate_grouped_data <- function(N = 2, bas = "piecewiseconstant", K = 50, beta = 0.2, sdtheta = 0.8){
  n <- N*K ### Total samples
  if(bas == "piecewiseconstant") {
    timelim <- baseline$timelim[1]
    tdom <- baseline$time
    haz <- baseline$hazard
    true <- data.frame(time = tdom, hazard = haz)
    u <- rnorm(K, sd = sdtheta)
    u <- rep(u, each = N)
    x <- rnorm(n, sd = 1)
    eta <- u + beta*x
    failtimes <- c()
    r <- runif(n)
    for (i in 1:n) {
      hazz <- haz * exp(eta[i])
      cumhaz <- cumsum(hazz*0.001)
      Surv <- exp(-cumhaz)
      Surv[1] <- 1
      failtimes[i] <- tdom[colSums(outer(Surv, r[i], `>`))]
    }
    data <- data_frame(x = x,times = failtimes, entry = rep(0,length(length(u))),censoring = ifelse(failtimes >= timelim,yes = 0, no=1))
    for (i in 1:length(data$censoring)) {
      if (data$censoring[i] == 1) {data$censoring[i] <- rbinom(n = 1,size = 1,prob = 0.9)}
    }
    data$group <- rep(1:K, each = N)
    data$true <- u
  }
  else if (bas == "regular") {
    timelim <- 200
    tdom <- seq(0, timelim, by = 0.001)
    haz <- rep(0, length(tdom))
    haz <- 0.2 * cos(0.15*tdom) + 0.3
    true <- data.frame(time = tdom, hazard = haz)
    u <- rnorm(K, sd = sdtheta)
    u <- rep(u, each = N)
    x <- rnorm(n, sd = 3)
    eta <- u + beta*x
    failtimes <- c()
    r <- runif(n)
    for (i in 1:n) {
      hazz <- haz * exp(eta[i])
      cumhaz <- cumsum(hazz*0.001)
      Surv <- exp(-cumhaz)
      Surv[1] <- 1
      failtimes[i] <- tdom[colSums(outer(Surv, r[i], `>`))]
    }
    data <- data_frame(x = x,times = failtimes, entry = rep(0,length(length(u))),censoring = ifelse(failtimes >= timelim,yes = 0, no=1))
    for (i in 1:length(data$censoring)) {
      if (data$censoring[i] == 1) {data$censoring[i] <- rbinom(n = 1,size = 1,prob = 0.9)}
    }
    data$group <- rep(1:K, each = N)
    data$true <- u
  }
  else if (bas == "constant") {
    timelim <- 200
    tdom <- seq(0, timelim, by = 0.001)
    haz <- rep(0.05, length(tdom))
    true <- data.frame(time = tdom, hazard = haz)
    u <- rnorm(K, sd = sdtheta)
    u <- rep(u, each = N)
    x <- rnorm(n, sd = 3)
    eta <- u + beta*x
    failtimes <- c()
    r <- runif(n)
    for (i in 1:n) {
      hazz <- haz * exp(eta[i])
      cumhaz <- cumsum(hazz*0.001)
      Surv <- exp(-cumhaz)
      Surv[1] <- 1
      failtimes[i] <- tdom[colSums(outer(Surv, r[i], `>`))]
    }
    data <- data_frame(x = x,times = failtimes, entry = rep(0,length(length(u))),censoring = ifelse(failtimes >= timelim,yes = 0, no=1))
    for (i in 1:length(data$censoring)) {
      if (data$censoring[i] == 1) {data$censoring[i] <- rbinom(n = 1,size = 1,prob = 0.9)}
    }
    data$group <- rep(1:K, each = N)
    data$true <- u
  }
  else{
    timelim <- 200
    tdom <- seq(0, timelim, by = 0.001)
    haz <- rep(0, length(tdom))
    cut <- 40
    for (i in 1:cut) {
      low <- as.numeric(quantile(tdom,(i-1)/cut))
      high <- as.numeric(quantile(tdom,(i)/cut))
      if(i %% 2 == 1){
          haz[tdom<=high & tdom > low] <- 0.01
      }
      else if(i %% 2 == 0){
        haz[tdom<=high & tdom > low] <- 0.25
      }
    }
    true <- data.frame(time = tdom, hazard = haz)
    u <- rnorm(K, sd = sdtheta)
    u <- rep(u, each = N)
    x <- rnorm(n, sd = 3)
    eta <- u + beta*x
    failtimes <- c()
    r <- runif(n)
    for (i in 1:n) {
      hazz <- haz * exp(eta[i])
      cumhaz <- cumsum(hazz*0.001)
      Surv <- exp(-cumhaz)
      Surv[1] <- 1
      failtimes[i] <- tdom[colSums(outer(Surv, r[i], `>`))] 
    }
    data <- data_frame(x = x,times = failtimes, entry = rep(0,length(length(u))),censoring = ifelse(failtimes >= timelim,yes = 0, no=1))
    for (i in 1:length(data$censoring)) {
      if (data$censoring[i] == 1) {data$censoring[i] <- rbinom(n = 1,size = 1,prob = 0.9)}
    }
    data$group <- rep(1:K, each = N)
    data$true <- u
  }
  data
}

prior.prec <- list(prec = list(prior = "pc.prec",
                               param = c(1, 0.5)))



# TMB function template
compile("03_coxph_frailty.cpp")
dyn.load(dynlib("03_coxph_frailty"))



###### Fit the proposed method
set.seed(1)
n <- K*N
data <- Simulate_grouped_data(N = N, bas = "piecewiseconstant", K = K, beta = beta, sdtheta = sd)
data <- abcoxph:::arrange_data(data)
dat <- tibble(x = data$x, t = data$times, cens = data$censoring, group = data$group)
dat$ranks <- rank(dat$t, ties.method = "min")
X <- as(as.matrix(dat$x),"dgTMatrix")
B <- as(abcoxph:::create_blist_element(u = dat$group)$B,"dgTMatrix")
D <- as(abcoxph:::create_diff_matrix(n), "dgTMatrix") ### n = K * N
### Setup TMB:
tmbdat <- list(
  # Design matrix (random and fixed)
  B = as(B,"dgTMatrix"),
  X = as(X,"dgTMatrix"),
  # Differencing matrix
  D = as(D,"dgTMatrix"),
  # Response
  ranks = as.integer(dat$ranks),
  cens = as.integer(dat$cens),
  # Prior params
  u = 1,
  alpha = 0.5,
  betaprec = 0.001)

tmbparams <- list(
  W = rep(0,ncol(B)+ncol(X)),
  theta = 0 # -2log(sigma)
)
##### Fitting:
ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "03_coxph_frailty",
  silent = TRUE
)
# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
# AGHQ
start_time <- Sys.time()
quad <- aghq::marginal_laplace_tmb(ff,15,0)
samps <- sample_marginal(quad,10000,interpolation = 'spline')
end_time <- Sys.time()
runtime_aghq <- end_time - start_time
#### runtime of the proposed method took 25 secs, each sample took 0.0025 secs


# Plot of theta posterior
prec_marg <- quad$marginals[[1]]
logpostsigma <- compute_pdf_and_cdf(prec_marg,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
with(logpostsigma,plot(transparam,pdf_transparam,type='l'))

theta_logprior <- function(theta, prior_alpha = tmbdat$alpha, 
                           prior_u = tmbdat$u) {
  lambda <- -log(prior_alpha)/prior_u
  log(lambda/2) - lambda * exp(-theta/2) - theta/2
}
priorfuncsigma <- function(x) (2/x) * exp(theta_logprior(-2 * log(x)))
aghq_hyper <- logpostsigma %>% select(c("transparam","pdf_transparam")) %>% mutate(method = "Proposed")
names(aghq_hyper)[1:2] <- c("x","y")
prior <- tibble(x = aghq_hyper$x, y = priorfuncsigma(aghq_hyper$x), method = "Prior")
hyper <- rbind(prior,aghq_hyper)
hyper %>% ggplot(aes(x,y)) +  geom_line(aes(linetype = method)) + xlim(c(0,2)) + xlab("SD") + ylab("density") + theme_classic(base_size = TEXT_SIZE)



#### Fit and Compare with STAN:
### STAN:
start_time <- Sys.time()
stanmod <- tmbstan(
  ff,
  chains = 4,
  cores = 4,
  iter = 35000,
  warmup = 25000,
  init = quad$optresults$mode,
  # init = 0,
  seed = 12345
)
end_time <- Sys.time()
runtime_MCMC <- end_time - start_time
#### Runtime for MCMC is 40.25 minutes; each iteration took 0.07 secs.
#### If run under the same time constraint as the proposed method, can only take 358 iterations.
summ <- summary(stanmod)$summary
beta_STAN <- summ[(K+1):(nrow(summ)-2),]
STAN_samples <- extract(stanmod)
theta_sample <- STAN_samples$theta
sd_sample <- sqrt(1/exp(theta_sample))

################### First case: look at the variance parameter for frailty

### KS statistics:
set.seed(100)
stansamps <- as.data.frame(stanmod)
numsamp <- nrow(stansamps)
quadsamp <- sample_marginal(quad,numsamp,interpolation = 'spline')$thetasamples[[1]]
stansamps$sigma <- exp(-stansamps$theta/2)
ks.test(stansamps$theta,quadsamp)$statistic
ks.test(stansamps$sigma,exp(-quadsamp/2))$statistic


# Look at the KS
# The distributions look pretty close:
hist(stansamps$theta,breaks = 100,freq=FALSE)
with(logpostsigma,lines(theta,pdf))
hist(stansamps$sigma,breaks = 100,freq=FALSE, xlim = c(0,1.5))
with(logpostsigma,lines(transparam,pdf_transparam))


############## Second Case: look at the single fixed parameter beta

post_means <- apply(samps$samps, 1, mean)
post_sds <- apply(samps$samps, 1, sd)
post_up <- apply(samps$samps, 1, quantile, probs = 0.975)
post_lo <- apply(samps$samps, 1, quantile, probs = 0.025)
post_sum_aghq <- data.frame(rbind(post_means,post_sds,post_up,post_lo))
fixed_effect <- tibble(AGHQ_mean = t(post_sum_aghq)[61,c(1)], 
                       AGHQ_sd = t(post_sum_aghq)[61,c(2)], 
                       STAN_mean = beta_STAN[1], STAN_sd = beta_STAN[3])


### compute KS for fixed effect:
beta_aghq <- samps$samps[61,]
beta_mcmc <- STAN_samples$W[,61]
ks.test(beta_mcmc,beta_aghq)$statistic

hist(beta_mcmc,breaks = 100,freq=FALSE)
with(density(beta_aghq, bw = 0.01),lines(x,y))
hist(stansamps$sigma,breaks = 100,freq=FALSE, xlim = c(0,1.5))
with(logpostsigma,lines(transparam,pdf_transparam))




############## Final Case: look at these 60 frailties, mean KS and max KS
KS_vec <- c()
for(i in 1:60){
  xi_aghq <- samps$samps[i,]
  xi_mcmc <- STAN_samples$W[,i]
  KS_vec[i] <- ks.test(xi_mcmc,xi_aghq)$statistic
}
mean(KS_vec)
max(KS_vec)
