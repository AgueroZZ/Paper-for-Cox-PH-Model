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

precompile()
TEXT_SIZE = 25

registerDoMC(ncores)


## Simulating function:
### Simulation Example
Simulate_baseline <- function(timelim = 300, breaks = 0.001, cut = 40){
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
        haz[tdom<=high & tdom > low] <-(0.05) *(c-min(c))
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
M <- 300

## simulate data:
set.seed(1234)
baseline <- Simulate_baseline()
ggplot(baseline,aes(x = time, y = hazard)) + geom_line() + theme(text = element_text(size = TEXT_SIZE))
ggsave(filename = "base1.png")

prior.prec <- list(prec = list(prior = "pc.prec",
                               param = c(tmbdat$u, tmbdat$a)))


compile("06_coxph_bspline.cpp")
dyn.load(dynlib("06_coxph_bspline"))



#### Aggregation through 300 aggregations:
do_once <- function(seed, truth = "complicated", N = 1000, baseline){
    set.seed(seed)
    data <- Simulate_data_extreme(baseline = baseline, truth = truth, N = N)
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
    # AGHQ
    quad <- aghq::marginal_laplace_tmb(ff,7,0)
    samps <- sample_marginal(quad,2000)
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
    ### INLA:
    prior.prec <- list(prec = list(prior = "pc.prec",
                                   param = c(tmbdat$u, tmbdat$a)))
    formula <- inla.surv(times,censoring) ~ f(exposure_binned,model = 'rw2',constr = T, hyper = prior.prec)
    Inlaresult <- inla(formula = formula, data = data, family = "coxph")
    fhat <- Inlaresult$summary.random$exposure_binned$mean
    fup <- Inlaresult$summary.random$exposure_binned$`0.975quant`
    flo <- Inlaresult$summary.random$exposure_binned$`0.025quant`
    plotINLA <- data.frame(x = Inlaresult$summary.random$exposure_binned$ID, f = fhat, up = fup, lo = flo)
    plotINLA$true <- truefunc(plotINLA$x) - mean(truefunc(plotINLA$x))
    rmse_inla <- sqrt( mean( ((plotINLA$f) - (plotINLA$true))^2 ) )
    mse_inla <- mean( ((plotINLA$f) - (plotINLA$true))^2 )
    
    covr_inla <- mean((plotINLA$true) <= (plotINLA$up) & (plotINLA$true) >= (plotINLA$lo))
    ## mgcv:
    mgcvmod_bs <- mgcv::gam(
      t ~ s(x,bs='bs',m=c(p-1,m),k=length(splineknots)-p),
      data = dat,
      family = cox.ph()
    )
    mgcv_bs_pred <- predict(mgcvmod_bs,newdata = data.frame(x = plotx),se.fit = TRUE)
    mgcv_bs_pred$fit <- mgcv_bs_pred$fit - mean(mgcv_bs_pred$fit)
    rmse_mgcv <- sqrt( mean( (((mgcv_bs_pred$fit) - ((truefunc(plotx)-mean(truefunc(plotx)))))^2)) )
    mse_mgcv <- mean( (((mgcv_bs_pred$fit) - ((truefunc(plotx)-mean(truefunc(plotx)))))^2))
    covr_mgcv <- mean((((truefunc(plotx)-mean(truefunc(plotx)))) <= (mgcv_bs_pred$fit + 2*mgcv_bs_pred$se.fit) & ((truefunc(plotx)-mean(truefunc(plotx)))) >= (mgcv_bs_pred$fit - 2*mgcv_bs_pred$se.fit)))
    result <- data.frame(
      method = c('AGHQ','MGCV','INLA'),
      rmse = c(rmse_aghq,rmse_mgcv,rmse_inla),
      coverage = c(covr_aghq,covr_mgcv,covr_inla),
      mse = c(mse_aghq, mse_mgcv, mse_inla)
    )
    result
}


### First Case: where baseline hazard is very complicated
result <- foreach(i = 1:M,.combine = rbind, .packages = c('foreach', 'stats', 'INLA', 'aghq', 'abcoxph')) %dopar% do_once(seed = i, truth = "complicated", N = 1000 , baseline)
agg_result <- result %>% group_by(method) %>% summarise(rmse = mean(rmse), coverage = mean(coverage), mse = mean(mse))
agg_result
save(result, file = "result_300_aggregations_complicated.rda")
save(agg_result, file = "aggresult_300_aggregations_complicated.rda")


#### Second Case: where baseline hazard is medium
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
baseline2 <- Simulate_baseline2()
ggplot(baseline2,aes(x = time, y = hazard)) + geom_line() + theme(text = element_text(size = TEXT_SIZE))
ggsave(filename = "base2.png")
result2 <- foreach(i = 1:M,.combine = rbind, .packages = c('foreach', 'stats', 'INLA', 'aghq', 'abcoxph')) %dopar% do_once(seed = i, truth = "complicated", N = 1000 , baseline2)
agg_result2 <- result2 %>% group_by(method) %>% summarise(rmse = mean(rmse), coverage = mean(coverage), mse = mean(mse))
agg_result2
save(result2, file = "result_300_aggregations_stepwise.rda")
save(agg_result2, file = "aggresult_300_aggregations_stepwise.rda")



#### Last Case: where baseline hazard is simple
set.seed(1234)
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
baseline3 <- Simulate_baseline3()
ggplot(baseline3,aes(x = time, y = hazard)) + geom_line() + theme(text = element_text(size = TEXT_SIZE))
ggsave(filename = "base3.png")
result3 <- foreach(i = 1:M,.combine = rbind, .packages = c('foreach', 'stats', 'INLA', 'aghq', 'abcoxph')) %dopar% do_once(seed = i, truth = "complicated", N = 1000 , baseline3)
agg_result3 <- result3 %>% group_by(method) %>% summarise(rmse = mean(rmse), coverage = mean(coverage), mse = mean(mse))
agg_result3
save(result3, file = "result_300_aggregations_stepwise.rda")
save(agg_result3, file = "aggresult_300_aggregations_stepwise.rda")

