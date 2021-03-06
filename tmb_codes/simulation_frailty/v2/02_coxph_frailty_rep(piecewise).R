### Coxph regression with frailty ###
library(tidyverse)
library(aghq)
library(mgcv)
library(Matrix)
library(TMB)
library(INLA)
library(tmbstan)
library(foreach)
library(doMC)
library(parallel)


precompile()
TEXT_SIZE = 25
ncores = 4
registerDoMC(ncores)
set.seed(123)

### Simulating function:
K = 100
beta = 0.2
sd = 0.8
M <- 300



Simulate_grouped_data <- function(N = 5, bas = "piecewiseconstant", K = 10, beta = 0.2, sdtheta = 0.8){
  n <- N*K ### Total samples
  if(bas == "piecewiseconstant") {
    timelim <- 200
    tdom <- seq(0, timelim, by = 0.001)
    haz <- rep(0, length(tdom))
    haz[1:round(length(haz)/3)] <- 0.01
    haz[(round(length(haz)/3) + 1):(round(1.5*length(haz))/3)] <- 0.25
    haz[(round(1.5*length(haz)/3) + 1):length(haz)] <- 0.5
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
    haz[1:round(length(haz)/3)] <- 0.01*sin(0.85*tdom[1:round(length(haz)/3)]) + 0.015
    haz[(round(length(haz)/3) + 1):(round(1.5*length(haz))/3)] <- 0.2 * cos(0.25*tdom[(round(length(haz)/3) + 1):(round(1.5*length(haz))/3)]) + 0.3
    haz[(round(1.5*length(haz)/3) + 1):length(haz)] <- 0.02*sin(0.55*tdom[(round(1.5*length(haz)/3) + 1):length(haz)]) + 0.05
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
# compile("02_coxph_frailty.cpp")
dyn.load(dynlib("02_coxph_frailty"))















######### Speeding up the simulation function:
do_once <- function(seed,beta, N, K, sd, bas = "constant"){
  set.seed(seed + rpois(1,10))
  n <- K*N
  data <- Simulate_grouped_data(N = N, bas = bas, K = K, beta = beta, sdtheta = sd)
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
    DLL = "02_coxph_frailty",
    silent = TRUE
  )
  # Hessian not implemented for RE models
  ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
  # AGHQ
  quad <- aghq::marginal_laplace_tmb(ff,15,0)
  samps <- sample_marginal(quad,2000)
  beta_est <- samps$samps[(K+1),]
  post_means <- apply(samps$samps, 1, mean)
  post_sds <- apply(samps$samps, 1, sd)
  post_up <- apply(samps$samps, 1, quantile, probs = 0.975)
  post_lo <- apply(samps$samps, 1, quantile, probs = 0.025)
  post_sum_aghq <- data.frame(rbind(post_means,post_sds,post_up,post_lo))
  rownames(post_sum_aghq) <- c("mean", "sd", "upper", "lower")
  colnames(post_sum_aghq) <- c(c(1:K),"beta")
  beta_cov_aghq <- ifelse(beta <= post_sum_aghq[,K+1][3] & beta >= post_sum_aghq[,K+1][4], 1, 0)
  beta_mse_aghq <- (post_sum_aghq[,K+1][1] - beta)^2
  ### INLA:
  prior.prec <- list(prec = list(prior = "pc.prec",
                                 param = c(tmbdat$u, tmbdat$a)))
  formula <- inla.surv(t,cens)~ x + f(group, model = "iid", hyper = prior.prec)
  Inlaresult <- inla(formula = formula, control.fixed = list(prec = 0.001), data = dat, family = "coxph")
  beta_cov_inla <- ifelse(beta <= Inlaresult$summary.fixed[2,]$'0.975quant' & beta >= Inlaresult$summary.fixed[2,]$'0.025quant', 1, 0)
  beta_mse_inla <- (Inlaresult$summary.fixed[2,]$mean - beta)^2
  frailty <- data %>% select(c(group,true)) %>% arrange(group) %>% unique(by = group)
  frailty$AGHQ <- as.numeric(post_sum_aghq[1,][-(K+1)])
  frailty$INLA <- as.numeric(Inlaresult$summary.random$group$mean)
  ### Random effects: MSE
  frailty_mse_aghq <- mean((frailty$true - frailty$AGHQ)^2)
  frailty_mse_inla <- mean((frailty$true - frailty$INLA)^2)
  ### Random effects: coverage
  frailty$AGHQ_up <- as.numeric(post_sum_aghq[3,][-(K+1)])
  frailty$INLA_up <- as.numeric(Inlaresult$summary.random$group$`0.975quant`)
  frailty$AGHQ_lo <- as.numeric(post_sum_aghq[4,][-(K+1)])
  frailty$INLA_lo <- as.numeric(Inlaresult$summary.random$group$`0.025quant`)
  ### AGHQ:
  frailty_cov_aghq <- mean((frailty$true >= frailty$AGHQ_lo & frailty$true <= frailty$AGHQ_up))
  ### INLA:
  frailty_cov_inla <- mean((frailty$true >= frailty$INLA_lo & frailty$true <= frailty$INLA_up))
  result <- c(beta_cov_aghq,beta_cov_inla,beta_mse_aghq,beta_mse_inla,frailty_cov_aghq,frailty_cov_inla,frailty_mse_aghq,frailty_mse_inla)
  result
}





### N = 2
time_begin <- Sys.time()
result2 <- foreach(i = 1:M,.combine = rbind, .packages = c('foreach', 'stats', 'INLA', 'aghq', 'abcoxph')) %dopar% do_once(seed = i, beta = beta, N = 2, K = K, sd = sd, bas = "piecewiseconstant")
time_end <- Sys.time()
time_end - time_begin
agg_means2 <- apply(result2, 2, mean)


## N = 4
time_begin <- Sys.time()
result4 <- foreach(i = 1:M,.combine = rbind, .packages = c('foreach', 'stats', 'INLA', 'aghq', 'abcoxph')) %dopar% do_once(seed = i, beta = beta, N = 4, K = K, sd = sd, bas = "piecewiseconstant")
time_end <- Sys.time()
time_end - time_begin
agg_means4 <- apply(result4, 2, mean)


### N = 6
time_begin <- Sys.time()
result6 <- foreach(i = 1:M,.combine = rbind, .packages = c('foreach', 'stats', 'INLA', 'aghq', 'abcoxph')) %dopar% do_once(seed = i, beta = beta, N = 6, K = K, sd = sd, bas = "piecewiseconstant")
time_end <- Sys.time()
time_end - time_begin
agg_means6 <- apply(result6, 2, mean)



## N = 8
time_begin <- Sys.time()
result8 <- foreach(i = 1:M,.combine = rbind, .packages = c('foreach', 'stats', 'INLA', 'aghq', 'abcoxph')) %dopar% do_once(seed = i, beta = beta, N = 8, K = K, sd = sd, bas = "piecewiseconstant")
time_end <- Sys.time()
time_end - time_begin
agg_means8 <- apply(result8, 2, mean)


## N = 10
time_begin <- Sys.time()
result10 <- foreach(i = 1:M,.combine = rbind, .packages = c('foreach', 'stats', 'INLA', 'aghq', 'abcoxph')) %dopar% do_once(seed = i, beta = beta, N = 10, K = K, sd = sd, bas = "piecewiseconstant")
time_end <- Sys.time()
time_end - time_begin
agg_means10 <- apply(result10, 2, mean)


### Combine:
aggresult <- rbind(agg_means2,agg_means4,agg_means6, agg_means8, agg_means10)
colnames(aggresult) <- c("beta_cov_aghq","beta_cov_inla","beta_mse_aghq","beta_mse_inla","frailty_cov_aghq","frailty_cov_inla", "frailty_mse_aghq", "frailty_mse_inla")
save(aggresult, file = "aggresultPiece.Rda")




