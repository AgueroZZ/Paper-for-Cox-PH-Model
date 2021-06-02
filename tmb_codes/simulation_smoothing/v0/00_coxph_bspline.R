### Coxph B-Spline regression ###
library(tidyverse)
library(aghq)
library(mgcv)
library(Matrix)
library(TMB)
library(INLA)
library(tmbstan)


## simulate data:
source("~/Documents/Paper-for-Cox-PH-Model/rcode for implementation/required_function.R")
set.seed(1234)
data <- Simulate_data(bas = "complicated", truth = "complicated", N = 500)
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
BD <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,m),k=length(splineknots)-p),data = dat,knots = list(x = splineknots))
BB <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,0),k=length(splineknots)-p),data = dat,knots = list(x = splineknots))
P <- BD$S[[1]] + BB$S[[1]] # O'sullivan spline "of the third kind"
P <- as(P,"dgTMatrix")
B <- BB$X
B <- as(B,"dgTMatrix")
X <- sparse.model.matrix(~1,data = dat)
D <- abcoxph:::create_diff_matrix(n)

### Setup TMB:
tmbdat <- list(
  # Design matrix
  BX = as(cbind(B,X),'dgTMatrix'),
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
  alpha = 0.5,
  # Beta precision
  betaprec = 0.001
)

tmbparams <- list(
  W = rep(0,ncol(B) + ncol(X)), # W = c(U,beta); U = B-Spline coefficients, beta = polynomial coefficients starting with intercept
  theta = 0 # -2log(sigma)
)




# TMB function template
compile("00_coxph_bspline.cpp")
dyn.load(dynlib("00_coxph_bspline"))


ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "00_coxph_bspline",
  silent = TRUE
)
# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)






# AGHQ
start_time <- Sys.time()
quad <- aghq::marginal_laplace_tmb(ff,7,0)
end_time <- Sys.time()
runtime_AGHQ <- end_time - start_time

# Plot of theta posterior
logpostsigma <- compute_pdf_and_cdf(quad$marginals[[1]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma,plot(transparam,pdf_transparam,type='l'))


# Inference for W
samps <- sample_marginal(quad,1000)

# Plot of curve
# Posterior mean
W <- apply(samps$samps,1,mean)
U <- W[1:ncol(P)]
beta <- W[(ncol(P)+1):length(W)]
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
samplestoplot <- samps$samps[ ,sample(1:ncol(samps$samps),100,replace = FALSE)]
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
samplestoplot <- samps_mcmc[sample(1:nrow(samps_mcmc),100,replace = FALSE),]
for (i in 1:ncol(samplestoplot)) {
  WS <- samplestoplot[i,]
  US <- WS[1:ncol(P)]
  plotetaS <- plotB %*% US
  plotetaS <- plotetaS - plotetaS[1]
  lines(plotx,plotetaS,col = 'lightgray')
}
lines(plotx,ploteta)
lines(plotx,truefunc(plotx)-truefunc(plotx)[1],col='red')






