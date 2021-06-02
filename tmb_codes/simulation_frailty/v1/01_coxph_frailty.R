### Coxph regression with frailty ###
library(tidyverse)
library(aghq)
library(mgcv)
library(Matrix)
library(TMB)
library(INLA)
library(tmbstan)
library(brinla)
library(ggplotify)
precompile()
TEXT_SIZE = 25


### Simulating function:
N = 2
K = 100
beta = 0.2
sd = 0.8
n = K*N
Simulate_grouped_data <- function(N = 5, bas = "constant", K = 10, beta = 0.2, sdtheta = 0.8){
  n <- N*K ### Total samples
  if(bas == "constant") {
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

## simulate data:
set.seed(123)
data <- Simulate_grouped_data(N = N, bas = "constant", K = K, beta = beta, sdtheta = sd)
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



# TMB function template
compile("01_coxph_frailty.cpp")
dyn.load(dynlib("01_coxph_frailty"))







##### Fitting:
ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "01_coxph_frailty",
  silent = TRUE
)
# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)

# AGHQ
quad <- aghq::marginal_laplace_tmb(ff,15,0)

# Plot of theta posterior
logpostsigma <- compute_pdf_and_cdf(quad$marginals[[1]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
with(logpostsigma,plot(transparam,pdf_transparam,type='l'))

# Inference for W
samps <- sample_marginal(quad,1000,interpolation = 'spline')
beta_est <- samps$samps[(K+1),]
hist(beta_est,breaks = 100)

### summary of estimates:
post_means <- apply(samps$samps, 1, mean)
post_sds <- apply(samps$samps, 1, sd)
post_up <- apply(samps$samps, 1, quantile, probs = 0.975)
post_lo <- apply(samps$samps, 1, quantile, probs = 0.025)
post_sum_aghq <- data.frame(rbind(post_means,post_sds,post_up,post_lo))
rownames(post_sum_aghq) <- c("mean", "sd", "upper", "lower")
colnames(post_sum_aghq) <- c(c(1:K),"beta")


### Using INLA:
prior.prec <- list(prec = list(prior = "pc.prec",
                               param = c(tmbdat$u, tmbdat$a)))
formula <- inla.surv(t,cens)~ x + f(group, model = "iid", hyper = prior.prec)
Inlaresult <- inla(formula = formula, control.fixed = list(prec = 0.001), data = dat, family = "coxph")
### beta:
Inlaresult$summary.fixed

### random effects:
Inlaresult$summary.random$group

### hyper-parameter:
brinla::bri.hyperpar.plot(Inlaresult, together = F)


####### Comparison:
### Fixed effect:
## Base R plot:
plot(x = Inlaresult$marginals.fixed$x[,1], y = Inlaresult$marginals.fixed$x[,2], col = "blue", type = 'l')
lines(x = seq(-0.1,1,by = 0.01), y = dnorm(seq(-0.1,1,by = 0.01), mean = post_sum_aghq[,(K+1)][1], sd = post_sum_aghq[,(K+1)][2]), col = "purple", type = "l", xlab = "beta", ylab = "density")
abline(v = beta, col = "red")
## ggplot:
fixplot <- data.frame(x = Inlaresult$marginals.fixed$x[,1], y = Inlaresult$marginals.fixed$x[,2], type = "INLA")
fixplot <- rbind(fixplot, data.frame(x = fixplot$x, y = dnorm(fixplot$x,mean = post_sum_aghq[,(K+1)][1], sd = post_sum_aghq[,(K+1)][2]), type = "AGHQ"))
fixplot %>% ggplot(aes(x = x, y = y, color = type)) + geom_line() + labs(x = "beta", y = "density", color = "Legend") +
  theme(text = element_text(size = TEXT_SIZE), legend.title = element_blank()) + geom_vline(xintercept = beta, color = "red") + 
  scale_color_manual(values=c("purple","blue"))
### Random effects:
frailty <- data %>% select(c(group,true)) %>% arrange(group) %>% unique(by = group)
frailty$AGHQ <- as.numeric(post_sum_aghq[1,][-(K+1)])
frailty$INLA <- as.numeric(Inlaresult$summary.random$group$mean)
frailty
### Random effects: MSE
mean((frailty$true - frailty$AGHQ)^2)
mean((frailty$true - frailty$INLA)^2)
### Plot: 95 credible intervals
frailty$AGHQ_up <- as.numeric(post_sum_aghq[3,][-(K+1)])
frailty$INLA_up <- as.numeric(Inlaresult$summary.random$group$`0.975quant`)
frailty$AGHQ_lo <- as.numeric(post_sum_aghq[4,][-(K+1)])
frailty$INLA_lo <- as.numeric(Inlaresult$summary.random$group$`0.025quant`)
ggplot(data = frailty, aes(x = group, y = true)) + geom_point() + geom_errorbar(aes(ymax = AGHQ_up, ymin = AGHQ_lo), color = "purple") + 
  geom_errorbar(aes(ymax = INLA_up, ymin = INLA_lo), color = "blue") + theme(text = element_text(size = TEXT_SIZE)) + ylab("Frailty")

#### Approximate coverage rate in each case:
### AGHQ:
mean((frailty$true >= frailty$AGHQ_lo & frailty$true <= frailty$AGHQ_up))
### INLA:
mean((frailty$true >= frailty$INLA_lo & frailty$true <= frailty$INLA_up))




#### Fit and Compare with STAN:
### STAN:
start_time <- Sys.time()
stanmod <- tmbstan(
  ff,
  chains = 4,
  cores = 4,
  iter = 2e04,
  warmup = 2000,
  # init = quad$optresults$mode,
  init = 0,
  seed = 123
)
end_time <- Sys.time()
runtime_MCMC <- end_time - start_time
summ <- summary(stanmod)$summary
beta_STAN <- summ[K+1,]
STAN_samples <- extract(stanmod)
theta_sample <- STAN_samples$theta
sd_sample <- sqrt(1/exp(theta_sample))

### Fixed effects comparison:
### Base R plot:
plot(x = Inlaresult$marginals.fixed$x[,1], y = Inlaresult$marginals.fixed$x[,2], col = "blue", type = 'l')
lines(x = seq(-0.1,1,by = 0.01), y = dnorm(seq(-0.1,1,by = 0.01), mean = post_sum_aghq[,(K+1)][1], sd = post_sum_aghq[,(K+1)][2]), col = "purple", type = "l", xlab = "beta", ylab = "density")
lines(x = seq(-0.1,1,by = 0.01), y = dnorm(seq(-0.1,1,by = 0.01), mean = summ[K+1,][1], sd = summ[K+1,][3]), col = "green", type = "l", xlab = "beta", ylab = "density")
abline(v = beta, col = "red")
### ggplot:
## ggplot:
fixplot <- data.frame(x = Inlaresult$marginals.fixed$x[,1], y = Inlaresult$marginals.fixed$x[,2], type = "INLA")
fixplot <- rbind(fixplot, data.frame(x = fixplot$x, y = dnorm(fixplot$x,mean = post_sum_aghq[,(K+1)][1], sd = post_sum_aghq[,(K+1)][2]), type = "AGHQ"))
fixplot <- rbind(fixplot, data.frame(x = fixplot$x, y = dnorm(fixplot$x,mean = summ[K+1,][1], sd = summ[K+1,][3]), type = "STAN"))
fixplot %>% ggplot(aes(x = x, y = y, color = type)) + geom_line() + labs(x = "beta", y = "density", color = "Legend") +
  theme(text = element_text(size = TEXT_SIZE), legend.title = element_blank()) + geom_vline(xintercept = beta, color = "red") +
  scale_color_manual(values=c("purple","blue",'green'))



### Hyperparameter comparison:
hyper_inla <- brinla::bri.hyperpar.plot(Inlaresult, together = F)
hyper_inla <- hyper_inla %>% filter(parameter == "SD for group")
hyper_stan <- density(sd_sample, bw = 0.1)
hyper_plot <- data.frame(x = logpostsigma$transparam, y = logpostsigma$pdf_transparam, type = "AGHQ")
hyper_plot <- rbind(hyper_plot, data.frame(x = hyper_inla[,1], y = hyper_inla[,2], type = "INLA"))
hyper_plot <- rbind(hyper_plot, data.frame(x = hyper_stan$x, y = hyper_stan$y, type = "STAN"))
hyper_plot %>% ggplot(aes(x,y,color = type)) + geom_line() + labs(x = "SD", y = "density", color = "Legend") +
  theme(text = element_text(size = TEXT_SIZE), legend.title = element_blank()) + geom_vline(xintercept = sd, color = "red") +
  scale_color_manual(values=c("purple","blue",'green')) + xlim(c(0,3))


### KS statistics:
stansamps <- as.data.frame(stanmod)
numsamp <- nrow(stansamps)
quadsamp <- sample_marginal(quad,numsamp,interpolation = 'spline')$thetasamples[[1]]
normsamp <- rnorm(numsamp,quad$optresults$mode,1/sqrt(quad$optresults$hessian))
ks.test(stansamps$theta,quadsamp)$statistic
ks.test(stansamps$theta,normsamp)$statistic

# Look at the KS
# The distributions look pretty close:
hist(stansamps$theta,breaks = 50,freq=FALSE)
with(logpostsigma,lines(theta,pdf))

# Compute the KS manually. Plot the ECDFs:
tt <- seq(-3,3,length.out=1e04)
quadecdf <- ecdf(quadsamp)(tt)
stanecdf <- ecdf(stansamps$theta)(tt)
plot(tt,quadecdf,type='l')
lines(tt,stanecdf,lty='dashed')

# KS is the max absolute difference:
theKS <- max(abs(stanecdf - quadecdf))
whereistheKS <- which.max(abs(stanecdf - quadecdf))
abline(v = tt[whereistheKS])
plot(tt,abs(stanecdf - quadecdf),type='l')


### Credible intervals for frailties
frailty$STAN <- STAN_samples$W[,1:K] %>% apply(MARGIN = 2, mean)
frailty$STAN_up <- STAN_samples$W[,1:K] %>% apply(MARGIN = 2, quantile, probs = 0.975)
frailty$STAN_lo <- STAN_samples$W[,1:K] %>% apply(MARGIN = 2, quantile, probs = 0.025)
colors <- c("INLA" = "blue", "AGHQ" = "purple", "STAN" = "green")
ggplot(data = frailty, aes(x = group, y = true)) + geom_point() + geom_errorbar(aes(ymax = AGHQ_up, ymin = AGHQ_lo, color = "AGHQ")) + 
  geom_errorbar(aes(ymax = INLA_up, ymin = INLA_lo, color = "INLA")) + 
  geom_errorbar(aes(ymax = STAN_up, ymin = STAN_lo, color = "STAN")) +
  labs(x = "group", y = "Frailty", color = "Legend") +
  theme(text = element_text(size = TEXT_SIZE), legend.title = element_blank()) +
  scale_color_manual(values = colors)
#### Approximate coverage rate in each case:
### AGHQ:
mean((frailty$true >= frailty$AGHQ_lo & frailty$true <= frailty$AGHQ_up))
### INLA:
mean((frailty$true >= frailty$INLA_lo & frailty$true <= frailty$INLA_up))
### STAN:
mean((frailty$true >= frailty$STAN_lo & frailty$true <= frailty$STAN_up))

### Random effects: MSE
mean((frailty$true - frailty$AGHQ)^2)
mean((frailty$true - frailty$INLA)^2)
mean((frailty$true - frailty$STAN)^2)










##### Repeat for 300 times:
set.seed(123)
repeating_frailties_M_times <- function(M = 100, beta, N, K, sd, bas = "constant"){
  result <- tibble()
  n <-N*K
  for(mm in 1:M){
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
      DLL = "01_coxph_frailty",
      silent = TRUE
    )
    # Hessian not implemented for RE models
    ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
    # AGHQ
    quad <- aghq::marginal_laplace_tmb(ff,3,0)
    samps <- sample_marginal(quad,1000)
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
    result <- rbind(result, c(beta_cov_aghq,beta_cov_inla,beta_mse_aghq,beta_mse_inla,frailty_cov_aghq,frailty_cov_inla,frailty_mse_aghq,frailty_mse_inla))
  }
  colnames(result) <- c("beta_cov_aghq","beta_cov_inla","beta_mse_aghq","beta_mse_inla","frailty_cov_aghq","frailty_cov_inla", "frailty_mse_aghq", "frailty_mse_inla")
  result
}
result <- repeating_frailties_M_times(M = 300, beta = beta, N = N, K = K, sd = sd, bas = "constant")
agg_means <- apply(result, 2, mean)


### Present a Table:
agg_table <- matrix(agg_means, ncol=2, nrow = 4, byrow = T)
colnames(agg_table) <- c("AGHQ","INLA")
rownames(agg_table) <- c("beta_cover","beta_mse", "frailty_cover", "frailty_mse")
agg_table



