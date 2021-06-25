### Coxph B-Spline regression ###
lib_loc <- '/home/ziang/lib'
library(tidyverse, lib = lib_loc)
library(aghq, lib = lib_loc)
library(mgcv)
library(Matrix)
library(rstan, lib = lib_loc)
library(TMB, lib = lib_loc)
library(INLA, lib = lib_loc)
library(tmbstan, lib = lib_loc)
library(foreach, lib = lib_loc)
library(doMC, lib = lib_loc)
library(parallel)
library(foreach)
library(abcoxph, lib = lib_loc)
library(mvQuad, lib = lib_loc)
library(survival)
library(brinla, lib = lib_loc)


precompile()
TEXT_SIZE = 25

ncores = 4
registerDoMC(ncores)
## Simulating function:
### Simulation Example
N = 1000
Simulate_data_extreme <- function(N = 1000, truth, RW2BINS = 50, baseline){
  tdom <- baseline$time
  timelim <- baseline$timelim[1]
  haz <- baseline$hazard
  if(truth == "smooth"){
    u <- runif(N)
    x <- runif(N,min = 0, max = 6)
    truefunc <- function(x) log((x + 1)^2) - 1
    eta <- truefunc(x)
  }
  else{
    u <- runif(N)
    x <- runif(N,min = -6, max = 6)
    truefunc <- function(x) 1.5*(sin(0.8*x))
    eta <- truefunc(x)
  }
  failtimes <- c()
  for (i in 1:N) {
    hazz <- haz * exp(eta[i])
    cumhaz <- cumsum(hazz*0.001)
    Surv <- exp(-cumhaz)
    failtimes[i] <- tdom[colSums(outer(Surv, u[i], `>`))]
  }
  data <- data_frame(x = x,times = failtimes, entry = rep(0,length(length(u))),censoring = ifelse(failtimes >= timelim,yes = 0, no=1))
  for (i in 1:length(data$censoring)) {
    if (data$censoring[i] == 1) {data$censoring[i] <- rbinom(n = 1,size = 1,prob = 0.9)}
  }
  data <- rename(data,exposure = x)
  data <- data %>% as_tibble() %>%
    mutate(exposure_binned = abcoxph:::bin_covariate(exposure,bins = RW2BINS,type = "equal"))
  data
}
construct_design <- function(x,splineknots,p,m) {
  BB <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,0),k=length(splineknots)-p),data = data.frame(x = x),knots = list(x = splineknots))
  B <- BB$X
  B <- as(B,"dgTMatrix")
  B
}
construct_penalty <- function(x,splineknots,p,m, noise = 0.0001) {
  BD <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,m),k=length(splineknots)-p),data = data.frame(x = x),knots = list(x = splineknots))
  BB <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,0),k=length(splineknots)-p),data = data.frame(x = x),knots = list(x = splineknots))
  # BD$S[[1]] + BB$S[[1]] # O'sullivan spline "of the third kind"
  BD$S[[1]] + diag(noise, ncol = ncol(BD$S[[1]]), nrow = nrow(BD$S[[1]])) #### add a very small noise to make the penalty matrix full rank
}
## simulate data:
prior.prec <- list(prec = list(prior = "pc.prec",
                               param = c(2, 0.5)))

### Compile C code:
compile("06_coxph_bspline.cpp")
dyn.load(dynlib("06_coxph_bspline"))

#### Proposed:
set.seed(1234)
Simulate_baseline2 <- function(timelim = 300, breaks = 0.001, cut = 30){
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
        haz[tdom<=high & tdom > low] <-0.05
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
baseline <- Simulate_baseline2()

set.seed(1)
data <- Simulate_data_extreme(baseline = baseline, truth = "complicated", N = N)
data <- abcoxph:::arrange_data(data)
dat <- tibble(x = data$exposure, t = data$times, cens = data$censoring)
dat$ranks <- rank(dat$t, ties.method = "min")
## setup smoothing part:
a <- min(dat$x)
b <- max(dat$x) # boundary
n <- nrow(dat)
# Order of spline
p <- 4 # 4 = cubic
# Order of derivative penalty
m <- 2
# Number of INTERIOR knots
d <- 46
# Number of knots
T <- d + p
# The knots
intknots <- seq(a,b,length.out = d)
leftknots <- seq(min(intknots)-(p-1),min(intknots)-1,by=1)
rightknots <- seq(max(intknots)+1,max(intknots)+p-1,by=1)
splineknots <- sort(unique(c(leftknots,intknots,rightknots)))
P <- as(construct_penalty(dat$x,splineknots,p,m,noise = 0.0001),'dgTMatrix')
B <- as(construct_design(dat$x,splineknots,p,m),'dgTMatrix')
D <- abcoxph:::create_diff_matrix(n)
### Setup TMB:
tmbdat <- list(
  # Design matrix
  BX = B,
  # Penalty matrix
  P = P,
  # Differencing matrix
  D = D,
  # Log determinant of penalty matrix (without the sigma part)
  logPdet = as.numeric(determinant(P,logarithm = TRUE)$modulus),
  # Response
  ranks = as.integer(dat$ranks),
  cens = as.integer(dat$cens),
  # Prior params
  u = 2,
  alpha = 0.5
)
tmbparams <- list(
  W = rep(0,ncol(B)), # W = c(U); U = B-Spline coefficients
  theta = 0 # -2log(sigma)
)
ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "06_coxph_bspline",
  silent = TRUE
)
# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
# Find proposed inference
start_time <- Sys.time()
quad <- aghq::marginal_laplace_tmb(ff,7,0)
samps <- sample_marginal(quad,2000)
end_time <- Sys.time()
runtime_aghq <- end_time - start_time
### Result from AGHQ:
W <- apply(samps$samps,1,mean)
U <- W[1:ncol(P)]
truefunc <- function(x) 1.5*(sin(0.8*x))
plotx <- seq(a,b,by=0.01)
sampled_lambda <- construct_design(plotx,splineknots,p,m) %*% samps$samps
for (i in 1:ncol(sampled_lambda)) {
  sampled_lambda[,i] <- sampled_lambda[,i] - mean(sampled_lambda[,i])
}
est_lambda <- apply(sampled_lambda,1,mean)
est_lambda <- est_lambda
lambda_lower <- apply(sampled_lambda,1,quantile,probs = .025)
lambda_upper <- apply(sampled_lambda,1,quantile,probs = .975)
true_lambda <- truefunc(plotx) - mean(truefunc(plotx))
rmse_aghq <- sqrt( mean( (est_lambda - true_lambda)^2 ) )
mse_aghq <- mean( (est_lambda - true_lambda)^2 )
covr_aghq <- mean(true_lambda <= lambda_upper & true_lambda >= lambda_lower)





################ MCMC:
#### Fit and Compare with STAN:
### STAN:
start_time <- Sys.time()
stanmod <- tmbstan(
  ff,
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 8000,
  init = quad$optresults$mode,
  # init = 0,
  seed = 12345
)
end_time <- Sys.time()
runtime_MCMC <- end_time - start_time

save(stanmod, file = "smoothMCMC10000.rda")
#### Runtime for 4000 iterations with 2000 warmups is 6.27997 hours
#### Now try for 10000 iterations with 8000 warmups


####### First case: KS statistics for the hyperparameter:
set.seed(100)
stansamps <- as.data.frame(stanmod)
numsamp <- nrow(stansamps)
quadsamp <- sample_marginal(quad,numsamp,interpolation = 'spline')$thetasamples[[1]]
stansamps$sigma <- exp(-stansamps$theta/2)
ks.test(stansamps$theta,quadsamp)$statistic
ks.test(stansamps$sigma,exp(-quadsamp/2))$statistic


############## Final Case: look at these , mean KS and max KS
STAN_samples <- extract(stanmod)
KS_vec <- c()
for(i in 1:ncol(P)){
  xi_aghq <- samps$samps[i,]
  xi_mcmc <- STAN_samples$W[,i]
  KS_vec[i] <- ks.test(xi_mcmc,xi_aghq)$statistic
}
mean(KS_vec)
max(KS_vec)

