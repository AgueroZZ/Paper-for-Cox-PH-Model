### Coxph B-Spline regression ###
library(tidyverse)
library(aghq)
library(mgcv)
library(Matrix)
library(TMB)
library(INLA)
library(tmbstan)
library(brinla)
precompile()

## Simulating function:
### Simulation Example
Simulate_data <- function(N = 500, bas, truth, RW2BINS = 50){
  if(bas == "constant") {
    timelim <- 200
    tdom <- seq(0, timelim, by = 0.001)
    haz <- rep(0, length(tdom))
    haz <- 0*tdom + 0.03
    true <- data.frame(time = tdom, hazard = haz)
    if(truth == "smooth") {
      u <- runif(N)
      x <- runif(N,min = 0, max = 6)
      truefunc <- function(x) log((x + 1)^2) - 1
      eta <- truefunc(x)
    }
    else{
      u <- runif(N)
      x <- runif(N,min = -5, max = 5)
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
  }
  else if (bas == "regular") {
    timelim <- 200
    tdom <- seq(0, timelim, by = 0.001)
    haz <- rep(0, length(tdom))
    haz <- 0.03 * cos(0.15*tdom) + 0.05
    true <- data.frame(time = tdom, hazard = haz)
    if(truth == "smooth"){
      u <- runif(N)
      x <- runif(N,min = 0, max = 6)
      truefunc <- function(x) log((x + 1)^2) - 1
      eta <- truefunc(x)
    }
    else{
      u <- runif(N)
      x <- runif(N,min = -5, max = 5)
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
  }
  else if (bas == "three_pieces") {
    timelim <- 300
    tdom <- seq(0, timelim, by = 0.001)
    haz <- rep(0, length(tdom))
    haz[1:round(length(haz)/3)] <- 0.01
    haz[(round(length(haz)/3) + 1):(round(1.5*length(haz))/3)] <- 0.25
    haz[(round(1.5*length(haz)/3) + 1):length(haz)] <- 0.5
    true <- data.frame(time = tdom, hazard = haz)
    if(truth == "smooth"){
      u <- runif(N)
      x <- runif(N,min = 0, max = 6)
      truefunc <- function(x) log((x + 1)^2) - 1
      eta <- truefunc(x)
    }
    else{
      u <- runif(N)
      x <- runif(N,min = -5, max = 5)
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
  }
  else if (bas == "extreme") {
    timelim <- 300
    tdom <- seq(0, timelim, by = 0.001)
    haz <- rep(0, length(tdom))
    cut <- 40
    for (i in 1:cut) {
      low <- as.numeric(quantile(tdom,(i-1)/cut))
      high <- as.numeric(quantile(tdom,(i)/cut))
      if(i %% 2 == 1){
        a <- runif(1,0,1)
        if(a > 0.3) haz[tdom<=high & tdom > low] <- 0.01
        else {
          c <- tdom[tdom<=high & tdom > low]
          haz[tdom<=high & tdom > low] <-(0.001) *(c-min(c))
        }
      }
      if(i %% 2 == 0){
        a <- runif(1,0,1)
        if(a > 0.8){
          c <- tdom[tdom<=high & tdom > low]
          haz[tdom<=high & tdom > low] <- 0.25
        }
        else{
          haz[tdom<=high & tdom > low] <- sample(c(0.001,0.002),size = 1,prob = c(0.5,0.5))
        }
      }
    }
    true <- data.frame(time = tdom, hazard = haz)
    if(truth == "smooth"){
      u <- runif(N)
      x <- runif(N,min = 0, max = 6)
      truefunc <- function(x) log((x + 1)^2) - 1
      eta <- truefunc(x)
    }
    else{
      u <- runif(N)
      x <- runif(N,min = -5, max = 5)
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
  }
  else{
    timelim <- 200
    tdom <- seq(0, timelim, by = 0.001)
    haz <- rep(0, length(tdom))
    haz[1:round(length(haz)/3)] <- 0.01*sin(0.85*tdom[1:round(length(haz)/3)]) + 0.02
    haz[(round(length(haz)/3) + 1):(round(2*length(haz))/3)] <- 0.03 * cos(0.25*tdom[(round(length(haz)/3) + 1):(round(2*length(haz))/3)]) + 0.05
    haz[(round(2*length(haz)/3) + 1):length(haz)] <- 0.01*sin(0.55*tdom[(round(2*length(haz)/3) + 1):length(haz)]) + 0.025
    true <- data.frame(time = tdom, hazard = haz)
    if(truth == "smooth"){
      u <- runif(N)
      x <- runif(N,min = 0, max = 6)
      truefunc <- function(x) log((x + 1)^2) - 1
      eta <- truefunc(x)
    }
    else{
      u <- runif(N)
      x <- runif(N,min = -5, max = 5)
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
  }
  data
}



## simulate data:
set.seed(1234)
data <- Simulate_data(bas = "extreme", truth = "complicated", N = 500)
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
d <- 42
# Number of knots
T <- d + 2*p
# The knots
intknots <- seq(a,b,length.out = d)
leftknots <- seq(min(intknots)-(p-1),min(intknots)-1,by=1)
rightknots <- seq(max(intknots)+1,max(intknots)+p-1,by=1)
splineknots <- sort(unique(c(leftknots,intknots,rightknots)))
construct_design <- function(x,splineknots,p,m) {
  BD <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,m),k=length(splineknots)-p),data = data.frame(x = x),knots = list(x = splineknots))
  BB <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,0),k=length(splineknots)-p),data = data.frame(x = x),knots = list(x = splineknots))
  P <- BD$S[[1]] + BB$S[[1]] # O'sullivan spline "of the third kind"
  P <- as(P,"dgTMatrix")
  B <- BB$X
  B <- as(B,"dgTMatrix")
  B
}
construct_penalty <- function(x,splineknots,p,m) {
  BD <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,m),k=length(splineknots)-p),data = data.frame(x = x),knots = list(x = splineknots))
  BB <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,0),k=length(splineknots)-p),data = data.frame(x = x),knots = list(x = splineknots))
  BD$S[[1]] + BB$S[[1]] # O'sullivan spline "of the third kind"
}

P <- as(construct_penalty(dat$x,splineknots,p,m),'dgTMatrix')
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
  u = 1,
  alpha = 0.5
)

tmbparams <- list(
  W = rep(0,ncol(B)), # W = c(U); U = B-Spline coefficients
  theta = 0 # -2log(sigma)
)

# TMB function template
compile("02_coxph_bspline.cpp")
dyn.load(dynlib("02_coxph_bspline"))


ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "02_coxph_bspline",
  silent = TRUE
)
# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)

# AGHQ
quad <- aghq::marginal_laplace_tmb(ff,7,0)

# Plot of theta posterior
logpostsigma <- compute_pdf_and_cdf(quad$marginals[[1]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma,plot(transparam,pdf_transparam,type='l'))

# Inference for W
samps <- sample_marginal(quad,1000)

# Plot of curve
# Posterior mean
W <- apply(samps$samps,1,mean)
U <- W[1:ncol(P)]
truefunc <- function(x) 1.5*(sin(0.8*x))

# Construct a plot
plotx <- seq(a,b,by=0.01)
plotdat <- data.frame(x = plotx)
plotB <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,0),k=length(splineknots)-p),data = plotdat,knots = list(x = splineknots))$X
# Construct the plot
ploteta <- plotB %*% U
ploteta <- ploteta - ploteta[1]
plot(plotx,ploteta,type='l',ylim=c(-4,4))
# Plot the samples
samplestoplot <- samps$samps[ ,sample(1:ncol(samps$samps),200,replace = FALSE)]
for (i in 1:ncol(samplestoplot)) {
  WS <- samplestoplot[ ,i]
  US <- WS[1:ncol(P)]
  plotetaS <- plotB %*% US
  plotetaS <- plotetaS - plotetaS[1]
  lines(plotx,plotetaS,col = 'lightgray')
}
lines(plotx,ploteta)
lines(plotx,truefunc(plotx)-truefunc(plotx)[1],col='red')


# get pointwise SD of eta
etaplot <- list()
for (i in 1:ncol(samps$samps)) {
  W <- samps$samps[ ,i]
  U <- W[1:ncol(P)]
  eta <- as.numeric(plotB %*% U)
  etaplot[[i]] <- data.frame(
    x = plotx,
    eta = eta - eta[1]
  )
}
etaplotframe <- purrr::reduce(etaplot,rbind) %>%
  group_by(x) %>%
  summarize(etamean = mean(eta),etasd = sd(eta)) %>%
  mutate(lower = etamean - 2*etasd,upper = etamean + 2*etasd)
with(etaplotframe,plot(x,etamean,type='l',ylim = c(-4,4)))
with(etaplotframe,lines(x,lower,type='l',lty='dashed'))
with(etaplotframe,lines(x,upper,type='l',lty='dashed'))
lines(plotx,truefunc(plotx)-truefunc(plotx)[1],col='red')

# Check mgcv
mgcvmod_bs <- mgcv::gam(
  t ~ s(x,bs='bs',m=c(p-1,m),k=length(splineknots)-p, pc = a),
  data = dat,
  family = cox.ph()
)
mgcv_bs_pred <- predict(mgcvmod_bs,newdata = data.frame(x = plotx),se.fit = TRUE)
plot(plotx,mgcv_bs_pred$fit,type = 'l',ylim = c(-4,4))
lines(plotx,(mgcv_bs_pred$fit - 2*mgcv_bs_pred$se.fit),type='l',lty='dashed')
lines(plotx,(mgcv_bs_pred$fit + 2*mgcv_bs_pred$se.fit),type='l',lty='dashed')
lines(plotx,truefunc(plotx)-truefunc(plotx)[1],col='red')


# Check INLA
cnsA1 <- matrix(rep(0,50),nrow = 1)
cnsA1[1] <- 1
conse <- matrix(0, nrow = 1, ncol = 1)

prior.prec <- list(prec = list(prior = "pc.prec",
                               param = c(tmbdat$u, tmbdat$a)))

formula <- inla.surv(times,censoring) ~ f(exposure_binned,model = 'rw2',constr = F, extraconstr = list(A=cnsA1,e=conse), hyper = prior.prec)
Inlaresult <- inla(formula = formula,control.inla = list(strategy = 'gaussian',int.strategy = 'grid'),data = data, family = "coxph")
fhat <- Inlaresult$summary.random$exposure_binned$mean
fup <- Inlaresult$summary.random$exposure_binned$`0.975quant`
flo <- Inlaresult$summary.random$exposure_binned$`0.025quant`
fhat[1] = 0
fup[1] = 0
flo[1] = 0
plotINLA <- data.frame(x = Inlaresult$summary.random$exposure_binned$ID, f = fhat, up = fup, lo = flo)
plotINLA$true <- truefunc(plotINLA$x) - truefunc(plotINLA$x)[1]
plot(plotINLA$x,plotINLA$f,type = 'l',ylim = c(-4,4))
lines(plotINLA$x,plotINLA$up,type='l',lty='dashed')
lines(plotINLA$x,plotINLA$lo,type='l',lty='dashed')
lines(plotINLA$x,plotINLA$true,col='red')





### STAN:
start_time <- Sys.time()
stanmod <- tmbstan(
  ff,
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 1000,
  init = quad$optresults$mode,
  seed = 123
)
end_time <- Sys.time()
runtime_MCMC <- end_time - start_time

summ <- summary(stanmod)
U_mcmc <- summ$summary[1:44,1]
ploteta_mcmc <- plotB %*% U_mcmc
ploteta_mcmc <- ploteta_mcmc - ploteta_mcmc[1]
plot(plotx,ploteta_mcmc,type='l',ylim=c(-4,4))


# Construct the plot
ploteta_mcmc <- plotB %*% U_mcmc
ploteta_mcmc <- ploteta_mcmc - ploteta_mcmc[1]
plot(plotx,ploteta_mcmc,type='l',ylim=c(-4,4))
# Plot the samples
samps_mcmc <- extract(stanmod)
samps_mcmc <- samps_mcmc$W
samplestoplot <- samps_mcmc[sample(1:nrow(samps_mcmc),1000,replace = FALSE),]
etaplot2 <- list()
for (i in 1:ncol(samplestoplot)) {
  WS <- samplestoplot[i,]
  US <- WS[1:ncol(P)]
  plotetaS <- plotB %*% US
  plotetaS <- plotetaS - plotetaS[1]
  lines(plotx,plotetaS,col = 'lightgray')
  etaplot2[[i]] <- data.frame(
    x = plotx,
    eta = plotetaS
  )
}
lines(plotx,ploteta)
lines(plotx,truefunc(plotx)-truefunc(plotx)[1],col='red')

etaplotframe2 <- purrr::reduce(etaplot2,rbind) %>%
  group_by(x) %>%
  summarize(etamean = mean(eta),etasd = sd(eta)) %>%
  mutate(lower = etamean - 2*etasd,upper = etamean + 2*etasd)
with(etaplotframe2,plot(x,etamean,type='l',ylim = c(-4,4)))
with(etaplotframe2,lines(x,lower,type='l',lty='dashed'))
with(etaplotframe2,lines(x,upper,type='l',lty='dashed'))
lines(plotx,truefunc(plotx)-truefunc(plotx)[1],col='red')



# RMSE and coverage. There is no "true" W; they are weights for the basis function
# representation. Evaluate the quality of the function approximation instead
sampled_lambda <- construct_design(plotx,splineknots,p,m) %*% samps$samps
for (i in 1:ncol(sampled_lambda)) {
  sampled_lambda[,i] <- sampled_lambda[,i] - sampled_lambda[1,i]
}
est_lambda <- apply(sampled_lambda,1,mean)
est_lambda <- est_lambda
lambda_lower <- apply(sampled_lambda,1,quantile,probs = .025)
lambda_upper <- apply(sampled_lambda,1,quantile,probs = .975)
true_lambda <- truefunc(plotx) - truefunc(plotx)[1]
plot(plotx,true_lambda,type='l',ylim = c(-5,5),col = 'red')
lines(plotx,est_lambda)
lines(plotx,lambda_lower,lty = 'dashed')
lines(plotx,lambda_upper,lty = 'dashed')
rmse_aghq <- sqrt( mean( (est_lambda - true_lambda)^2 ) )
covr_aghq <- mean(true_lambda <= lambda_upper & true_lambda >= lambda_lower)


### Check STAN:
stansamps <- as.data.frame(stanmod)
# sigma
hist(exp(-stansamps$theta/2),freq=FALSE,breaks = 100)
with(logpostsigma,lines(transparam,pdf_transparam)) # AGHQ
ss <- seq(1,6,length.out = 1e03)
tt <- -2*log(ss)
with(quad,lines(ss,dnorm(tt,optresults$mode,1/sqrt(optresults$hessian))*abs(numDeriv::grad(function(x) -2*log(x),ss)),lty='dashed'))
# ks disance
numsamp <- nrow(stansamps)
quadsamp <- sample_marginal(quad,numsamp)$thetasamples[[1]]
normsamp <- rnorm(numsamp,quad$optresults$mode,1/sqrt(quad$optresults$hessian))
ks.test(stansamps$theta,quadsamp)$statistic
ks.test(stansamps$theta,normsamp)$statistic




# metrics
plot(plotx,(mgcv_bs_pred$fit),type = 'l',ylim = c(-5,5))
lines(plotx,(mgcv_bs_pred$fit - 2*mgcv_bs_pred$se.fit),type='l',lty='dashed')
lines(plotx,(mgcv_bs_pred$fit + 2*mgcv_bs_pred$se.fit),type='l',lty='dashed')
lines(plotx,(truefunc(plotx)-truefunc(plotx)[1]),col='red')
rmse_mgcv <- sqrt( mean( (((mgcv_bs_pred$fit) - ((truefunc(plotx)-truefunc(plotx)[1])))^2)[-1]) )
covr_mgcv <- mean((((truefunc(plotx)-truefunc(plotx)[1])) <= (mgcv_bs_pred$fit + 2*mgcv_bs_pred$se.fit) & ((truefunc(plotx)-truefunc(plotx)[1])) >= (mgcv_bs_pred$fit - 2*mgcv_bs_pred$se.fit))[-1])


plot(plotINLA$x,(plotINLA$f),type = 'l',ylim = c(-5,5))
lines(plotINLA$x,(plotINLA$lo),type='l',lty='dashed')
lines(plotINLA$x,(plotINLA$up),type='l',lty='dashed')
lines(plotINLA$x,plotINLA$true,col='red')
rmse_inla <- sqrt( mean( ((plotINLA$f) - (plotINLA$true))^2 ) )
covr_inla <- mean((plotINLA$true) <= (plotINLA$up) & (plotINLA$true) >= (plotINLA$lo))
data.frame(
  method = c('AGHQ','MGCV','INLA'),
  rmse = c(rmse_aghq,rmse_mgcv,rmse_inla),
  coverage = c(covr_aghq,covr_mgcv,covr_inla)
)
















#### Aggregation through 300 aggregations:
set.seed(1234)
repeating_smoothing_M_times <- function(M = 100, truth = "complicated", N = 500 , bas = "extreme"){
  result <- tibble()
  for(mm in 1:M){
    data <- Simulate_data(bas = "extreme", truth = "complicated", N = 500)
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
    d <- 42
    # Number of knots
    T <- d + 2*p
    # The knots
    intknots <- seq(a,b,length.out = d)
    leftknots <- seq(min(intknots)-(p-1),min(intknots)-1,by=1)
    rightknots <- seq(max(intknots)+1,max(intknots)+p-1,by=1)
    splineknots <- sort(unique(c(leftknots,intknots,rightknots)))
    construct_design <- function(x,splineknots,p,m) {
      BD <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,m),k=length(splineknots)-p),data = data.frame(x = x),knots = list(x = splineknots))
      BB <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,0),k=length(splineknots)-p),data = data.frame(x = x),knots = list(x = splineknots))
      P <- BD$S[[1]] + BB$S[[1]] # O'sullivan spline "of the third kind"
      P <- as(P,"dgTMatrix")
      B <- BB$X
      B <- as(B,"dgTMatrix")
      B
    }
    construct_penalty <- function(x,splineknots,p,m) {
      BD <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,m),k=length(splineknots)-p),data = data.frame(x = x),knots = list(x = splineknots))
      BB <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,0),k=length(splineknots)-p),data = data.frame(x = x),knots = list(x = splineknots))
      BD$S[[1]] + BB$S[[1]] # O'sullivan spline "of the third kind"
    }
    
    P <- as(construct_penalty(dat$x,splineknots,p,m),'dgTMatrix')
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
      u = 1,
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
      DLL = "02_coxph_bspline",
      silent = TRUE
    )
    # Hessian not implemented for RE models
    ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
    # AGHQ
    quad <- aghq::marginal_laplace_tmb(ff,7,0)
    samps <- sample_marginal(quad,1000)
    W <- apply(samps$samps,1,mean)
    U <- W[1:ncol(P)]
    truefunc <- function(x) 1.5*(sin(0.8*x))
    plotx <- seq(a,b,by=0.01)
    sampled_lambda <- construct_design(plotx,splineknots,p,m) %*% samps$samps
    for (i in 1:ncol(sampled_lambda)) {
      sampled_lambda[,i] <- sampled_lambda[,i] - sampled_lambda[1,i]
    }
    est_lambda <- apply(sampled_lambda,1,mean)
    est_lambda <- est_lambda
    lambda_lower <- apply(sampled_lambda,1,quantile,probs = .025)
    lambda_upper <- apply(sampled_lambda,1,quantile,probs = .975)
    true_lambda <- truefunc(plotx) - truefunc(plotx)[1]
    rmse_aghq <- sqrt( mean( (est_lambda - true_lambda)^2 ) )
    covr_aghq <- mean(true_lambda <= lambda_upper & true_lambda >= lambda_lower)
    ### INLA:
    cnsA1 <- matrix(rep(0,50),nrow = 1)
    cnsA1[1] <- 1
    conse <- matrix(0, nrow = 1, ncol = 1)
    prior.prec <- list(prec = list(prior = "pc.prec",
                                   param = c(tmbdat$u, tmbdat$a)))
    formula <- inla.surv(times,censoring) ~ f(exposure_binned,model = 'rw2',constr = F, extraconstr = list(A=cnsA1,e=conse), hyper = prior.prec)
    Inlaresult <- inla(formula = formula,control.inla = list(strategy = 'gaussian',int.strategy = 'grid'),data = data, family = "coxph")
    fhat <- Inlaresult$summary.random$exposure_binned$mean
    fup <- Inlaresult$summary.random$exposure_binned$`0.975quant`
    flo <- Inlaresult$summary.random$exposure_binned$`0.025quant`
    fhat[1] = 0
    fup[1] = 0
    flo[1] = 0
    plotINLA <- data.frame(x = Inlaresult$summary.random$exposure_binned$ID, f = fhat, up = fup, lo = flo)
    plotINLA$true <- truefunc(plotINLA$x) - truefunc(plotINLA$x)[1]
    rmse_inla <- sqrt( mean( ((plotINLA$f) - (plotINLA$true))^2 ) )
    covr_inla <- mean((plotINLA$true) <= (plotINLA$up) & (plotINLA$true) >= (plotINLA$lo))
    ## mgcv:
    mgcvmod_bs <- mgcv::gam(
      t ~ s(x,bs='bs',m=c(p-1,m),k=length(splineknots)-p, pc = a),
      data = dat,
      family = cox.ph()
    )
    mgcv_bs_pred <- predict(mgcvmod_bs,newdata = data.frame(x = plotx),se.fit = TRUE)
    rmse_mgcv <- sqrt( mean( (((mgcv_bs_pred$fit) - ((truefunc(plotx)-truefunc(plotx)[1])))^2)[-1]) )
    covr_mgcv <- mean((((truefunc(plotx)-truefunc(plotx)[1])) <= (mgcv_bs_pred$fit + 2*mgcv_bs_pred$se.fit) & ((truefunc(plotx)-truefunc(plotx)[1])) >= (mgcv_bs_pred$fit - 2*mgcv_bs_pred$se.fit))[-1])
    resultnew <- data.frame(
      method = c('AGHQ','MGCV','INLA'),
      rmse = c(rmse_aghq,rmse_mgcv,rmse_inla),
      coverage = c(covr_aghq,covr_mgcv,covr_inla)
    )
    result <- rbind(result,resultnew)
  }
  result
}
result <- repeating_smoothing_M_times(M = 300, truth = "complicated", N = 500, bas = "extreme")
agg_result <- result %>% group_by(method) %>% summarise(rmse = mean(rmse), coverage = mean(coverage))
agg_result











