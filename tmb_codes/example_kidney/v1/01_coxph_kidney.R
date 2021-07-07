# ### CoxPH regression for the kidney example ###
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
# library(survival)
# library(brinla, lib = lib_loc)

### CoxPH regression for the kidney example ###

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
library(survival)
library(brinla)


precompile()
TEXT_SIZE = 25



## tidy up the data:
data <- survival::kidney
data <- data.frame(times = data$time, status = data$status, age = data$age,sex = data$sex,GN = ifelse(data$disease == "GN",1,0),AN = ifelse(data$disease == "AN",1,0),PKD = ifelse(data$disease == "PKD",1,0), id = data$id)
data <- abcoxph:::arrange_data(data)
data$ranks <- rank(data$times, ties.method = "min")
X <- as(as.matrix(data[,c(3:7)]),"dgTMatrix")
B <- as(abcoxph:::create_blist_element(u = data$id)$B,"dgTMatrix")
n <- nrow(data)
D <- as(abcoxph:::create_diff_matrix(n), "dgTMatrix") ### n = K * N


### Setup TMB:
tmbdat <- list(
  # Design matrix (random and fixed)
  B = as(B,"dgTMatrix"),
  X = as(X,"dgTMatrix"),
  # Differencing matrix
  D = as(D,"dgTMatrix"),
  # Response
  ranks = as.integer(data$ranks),
  cens = as.integer(data$status),
  # Prior params
  u = 2,
  alpha = 0.5,
  betaprec = 0.001)

tmbparams <- list(
  W = rep(0,ncol(B)+ncol(X)),
  theta = 0 # -2log(sigma)
)



# TMB function template
# compile("01_coxph_kidney.cpp")
dyn.load(dynlib("01_coxph_kidney"))







##### Fitting:
ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "01_coxph_kidney",
  silent = TRUE
)
# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)

# AGHQ
start_time <- Sys.time()
quad <- aghq::marginal_laplace_tmb(ff,15,0)

# Plot of theta posterior
prec_marg <- quad$marginals[[1]]
logpostsigma <- compute_pdf_and_cdf(prec_marg,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
with(logpostsigma,plot(transparam,pdf_transparam,type='l'))

# Inference for W
K <- ncol(B)
samps <- sample_marginal(quad,10000,interpolation = 'spline')
beta_est <- samps$samps[(K+1):nrow(samps$samps),]
end_time <- Sys.time()
runtime_aghq <- end_time - start_time
### summary of estimates:
post_means <- apply(samps$samps, 1, mean)
post_sds <- apply(samps$samps, 1, sd)
post_up <- apply(samps$samps, 1, quantile, probs = 0.975)
post_lo <- apply(samps$samps, 1, quantile, probs = 0.025)
post_sum_aghq <- data.frame(rbind(post_means,post_sds,post_up,post_lo))
rownames(post_sum_aghq) <- c("mean", "sd", "upper", "lower")
colnames(post_sum_aghq) <- c(c(1:K),"age","sex","GN","AN","PKD")

### Using INLA:
prior.prec <- list(prec = list(prior = "pc.prec",
                               param = c(tmbdat$u, tmbdat$a)))
formula <- inla.surv(times,status)~ age + sex + GN + AN + PKD + f(id, model = "iid", hyper = prior.prec)
Inlaresult <- inla(formula = formula, control.fixed = list(prec = 0.001), data = data, family = "coxph", control.compute=list(config = TRUE))
### beta:
Inlaresult$summary.fixed

### random effects:
Inlaresult$summary.random$id

### hyper-parameter:
inla_hyper <- brinla::bri.hyperpar.plot(Inlaresult, together = F)
inla_hyper <- inla_hyper %>% filter(parameter == "SD for id") %>% select(1:2) %>% mutate(method = "INLA")
aghq_hyper <- logpostsigma %>% select(c("transparam","pdf_transparam")) %>% mutate(method = "Proposed")
names(aghq_hyper)[1:2] <- c("x","y")
hyper <- rbind(inla_hyper,aghq_hyper)
hyper %>% ggplot(aes(x,y,type = method)) + geom_line() + xlim(c(0,2)) + xlab("SD") + ylab("density") + theme_classic(base_size = TEXT_SIZE)


theta_logprior <- function(theta, prior_alpha = tmbdat$alpha, 
                                      prior_u = tmbdat$u) {
  lambda <- -log(prior_alpha)/prior_u
  log(lambda/2) - lambda * exp(-theta/2) - theta/2
}
priorfuncsigma <- function(x) (2/x) * exp(theta_logprior(-2 * log(x)))
prior <- tibble(x = hyper$x, y = priorfuncsigma(hyper$x), method = "Prior")
hyper <- rbind(prior,hyper)
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
summ <- summary(stanmod)$summary
beta_STAN <- summ[(K+1):(nrow(summ)-2),]
STAN_samples <- extract(stanmod)
theta_sample <- STAN_samples$theta
sd_sample <- sqrt(1/exp(theta_sample))



### KS statistics:
stansamps <- as.data.frame(stanmod)
numsamp <- nrow(stansamps)
quadsamp <- sample_marginal(quad,numsamp,interpolation = 'spline')$thetasamples[[1]]
normsamp <- rnorm(numsamp,quad$optresults$mode,1/sqrt(quad$optresults$hessian))
stansamps$sigma <- exp(-stansamps$theta/2)
ks.test(stansamps$theta,quadsamp)$statistic
ks.test(stansamps$sigma,exp(-quadsamp/2))$statistic

# Look at the KS
# The distributions look pretty close:
hist(stansamps$theta,breaks = 50,freq=FALSE)
with(logpostsigma,lines(theta,pdf))
hist(stansamps$sigma,breaks = 50,freq=FALSE, xlim = c(0,3), ylim = c(0,1.5))
with(logpostsigma,lines(transparam,pdf_transparam))

# Compute the KS manually. Plot the ECDFs:
tt <- seq(-3,6,length.out=1e04)
quadecdf <- ecdf(quadsamp)(tt)
stanecdf <- ecdf(stansamps$theta)(tt)
plot(tt,quadecdf,type='l')
lines(tt,stanecdf,lty='dashed')

# KS is the max absolute difference:
theKS <- max(abs(stanecdf - quadecdf))
whereistheKS <- which.max(abs(stanecdf - quadecdf))
abline(v = tt[whereistheKS])
plot(tt,abs(stanecdf - quadecdf),type='l')


ggplot(stansamps, aes(x = sigma)) + 
  geom_histogram(fill = "gray", color = "gray", stat = "density") + 
  geom_line(data = logpostsigma, aes(x = transparam, y = pdf_transparam)) + xlim(0,3) +
  geom_line(data = inla_hyper, aes(x,y), linetype = "dashed") + xlab("SD") + ylab("Density") +
  theme_classic(base_size = TEXT_SIZE) 



### Compare fixed effect estimates:
fixed_effect <- tibble(AGHQ_mean = t(post_sum_aghq)[c(39:43),c(1)], 
                       AGHQ_sd = t(post_sum_aghq)[c(39:43),c(2)], 
                       STAN_mean = beta_STAN[,1], STAN_sd = beta_STAN[,3],
                       INLA_mean = Inlaresult$summary.fixed[-1,1], 
                       INLA_sd = Inlaresult$summary.fixed[-1,2])




############# Fixed effect KS:
KS_vec_fx <- c()
for(i in 39:43){
  fx_aghq <- samps$samps[i,]
  fx_mcmc <- STAN_samples$W[,i]
  KS_vec_fx <- c(KS_vec_fx,ks.test(fx_mcmc,fx_aghq)$statistic)
}
mean(KS_vec_fx)
max(KS_vec_fx)


########### Frailties KS:
KS_vec <- c()
for(i in 1:38){
  xi_aghq <- samps$samps[i,]
  xi_mcmc <- STAN_samples$W[,i]
  KS_vec[i] <- ks.test(xi_mcmc,xi_aghq)$statistic
}
mean(KS_vec)
max(KS_vec)




############## Sampling from INLA:
set.seed(123)
inlasample <- inla.posterior.sample(10000, result = Inlaresult)

### hyper
inla.varpara.sample <- numeric(10000)
for (i in 1:10000) {
  inla.varpara.sample[i] <- inlasample[[i]]$hyperpar[1]
}
ks.test(stansamps$sigma,1/sqrt(inla.varpara.sample))$statistic

### fixed:
inla.fx.sample <- matrix(nrow = 5, ncol = 10000)
for (i in 1:10000) {
  for (j in 1:5) {
    inla.fx.sample[j,i] <- inlasample[[i]]$latent[303 + j]
  }
}
KS_vec <- c()
for(i in 1:5){
  fixed_inla <- inla.fx.sample[i,]
  ii <- 38 + i
  fixed_mcmc <- STAN_samples$W[,ii]
  KS_vec[i] <- ks.test(fixed_mcmc,fixed_inla)$statistic
}
mean(KS_vec)
max(KS_vec)


########### Frailties KS:
inla.rd.sample <- matrix(nrow = 38, ncol = 10000)
for (i in 1:10000) {
  for (j in 1:38) {
    inla.rd.sample[j,i] <- inlasample[[i]]$latent[248 + j]
  }
}
KS_vec <- c()
for(i in 1:38){
  xi_inla <- inla.rd.sample[i,]
  xi_mcmc <- STAN_samples$W[,i]
  KS_vec[i] <- ks.test(xi_mcmc,xi_inla)$statistic
}
mean(KS_vec)
max(KS_vec)
