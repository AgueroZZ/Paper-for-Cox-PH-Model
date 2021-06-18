### Poisson B-Spline regression ###

library(tidyverse)
library(aghq)
library(mgcv)
library(Matrix)
library(TMB)
precompile()
library(INLA)

library(tmbstan)
library(parallel)
options(mc.cores = parallel::detectCores())

globalpath <- "/storage/phd/projects/random-walks/"
tmbpath <- paste0(globalpath,"code/tmb/")

set.seed(8079)
# with n = 1000 observed numerical instability in INLA
# did not observe with n = 100.
n <- 100 # sample size
a <- 0
b <- 1 # boundary

x <- runif(n,a,b) # covariate

truefunc <- function(z) log(1 + sin(2*pi*z^2)) # log scale
lam <- exp(truefunc(x))
y <- rpois(n,lam)
plot(seq(a,b,length.out = 1e03),exp(truefunc(seq(a,b,length.out = 1e03))),type='l',ylim = range(y))
# points(x,lam,pch=20)
points(x,y,pch=4)

dat <- data.frame(x = x,y = y)

## Setup TMB ----

# Order of spline
p <- 4 # 4 = cubic
# Order of derivative penalty
m <- 2
# Number of INTERIOR knots
d <- 12
# Number of knots
T <- d + 2*p
# The knots
intknots <- seq(a,b,length.out = d)
leftknots <- seq(min(intknots)-(p-1),min(intknots)-1,by=1)
rightknots <- seq(max(intknots)+1,max(intknots)+p-1,by=1)
splineknots <- sort(unique(c(leftknots,intknots,rightknots)))
plot(splineknots,exp(truefunc(splineknots)),pch = 20,ylim = range(y))
lines(seq(min(splineknots),max(splineknots),by=.01),exp(truefunc(seq(min(splineknots),max(splineknots),by=.01))))
points(x,y,pch = 4)

# # Define the spline basis
# splineBasis <- fda::create.bspline.basis(rangeval = range(splineknots),norder = p,breaks = splineknots)
# # Design matrix
# B <- fda::bsplineS(x = x,breaks = splineknots,norder = p,returnMatrix = TRUE)
# # Penalty matrix
# P <- fda::bsplinepen(splineBasis,Lfdobj = m)


# Or, with mgcv?
construct_design <- function(x,splineknots,p,m) {
  BD <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,m),k=length(splineknots)-p),data = data.frame(x = x),knots = list(x = splineknots))
  BB <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,0),k=length(splineknots)-p),data = data.frame(x = x),knots = list(x = splineknots))
  
  # rankMatrix(BD$S[[1]]) # 10
  # rankMatrix(BB$S[[1]]) # 12
  # rankMatrix(BD$S[[1]] + BB$S[[1]]) # 12
  
  P <- BD$S[[1]] + BB$S[[1]] # O'sullivan spline "of the third kind"
  P <- as(P,"dgTMatrix")
  B <- BB$X
  B <- as(B,"dgTMatrix")
  
  
  # Polynomial basis
  # degree m-1 (order m) orthogonal polynomial
  # X <- sparse.model.matrix(~poly(x,degree = m-1,raw = FALSE))
  X <- sparse.model.matrix(~1,data = data.frame(x = x))
  cbind(B,X)
}

construct_penalty <- function(x,splineknots,p,m) {
  BD <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,m),k=length(splineknots)-p),data = data.frame(x = x),knots = list(x = splineknots))
  BB <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,0),k=length(splineknots)-p),data = data.frame(x = x),knots = list(x = splineknots))
  
  BD$S[[1]] + BB$S[[1]] # O'sullivan spline "of the third kind"
}

P <- as(construct_penalty(dat$x,splineknots,p,m),'dgTMatrix')
design <- as(construct_design(dat$x,splineknots,p,m),'dgTMatrix')

tmbdat <- list(
  # Design matrix
  BX = design,
  # Penalty matrix
  P = P,
  # Log determinant of penalty matrix (without the sigma part)
  logPdet = as.numeric(determinant(P,logarithm = TRUE)$modulus),
  # Response
  y = y,
  # Prior params
  u = 1,
  alpha = 0.5,
  # Beta precision
  betaprec = 0.001
)

tmbparams <- list(
  W = rep(0,ncol(design)), # W = c(U,beta); U = B-Spline coefficients, beta = polynomial coefficients starting with intercept
  theta = 0 # -2log(sigma)
)

# TMB function template
# compile(paste0(tmbpath,"03_poisson_bspline.cpp"))
dyn.load(dynlib(paste0(globalpath,"code/tmb/03_poisson_bspline")))

ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "03_poisson_bspline",
  silent = TRUE
)
# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)

system.time(ff$fn(1))
system.time(ff$gr(1))
system.time(ff$he(1))


# AGHQ
quad <- aghq::marginal_laplace_tmb(ff,7,0)


# Plot of theta posterior
logpostsigma <- compute_pdf_and_cdf(quad$marginals[[1]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)))
with(logpostsigma,plot(transparam,pdf_transparam,type='l'))

# Inference for W
samps <- sample_marginal(quad,1e03)



# Plot of curve
# Posterior mean
W <- apply(samps$samps,1,mean)
U <- W[1:ncol(P)]
beta <- W[(ncol(P)+1):length(W)]

# Construct a plot
plotx <- seq(a,b,by=0.01)
plotdat <- data.frame(x = plotx)
plotB <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,0),k=length(splineknots)-p),data = plotdat,knots = list(x = splineknots))$X
# plotX <- sparse.model.matrix(~poly(plotx,degree = m-1,raw = FALSE))
plotX <- sparse.model.matrix(~1,data = plotdat)
# Construct the plot
ploteta <- plotB %*% U + plotX %*% beta
plot(plotx,exp(ploteta),type='l',ylim=c(0,12))
# Plot the samples
samplestoplot <- samps$samps[ ,sample(1:ncol(samps$samps),100,replace = FALSE)]
for (i in 1:ncol(samplestoplot)) {
  WS <- samplestoplot[ ,i]
  US <- WS[1:ncol(P)]
  betaS <- WS[(ncol(P)+1):length(WS)]
  plotetaS <- plotB %*% US + plotX %*% betaS
  lines(plotx,exp(plotetaS),col = 'lightgray')
}
lines(plotx,exp(ploteta))
lines(plotx,exp(truefunc(plotx)),col='red')
points(x,y,pch=4)

# NOTE: for some reason, the underlying polynomial BLOWS IT UP
# HAVE TO CHECK THIS
# It did NOT work with linear
# It DID work with intercept only.
# Note with intercept, still have sum(U) =/= 0

# RMSE and coverage. There is no "true" W; they are weights for the basis function
# representation. Evaluate the quality of the function approximation instead
sampled_lambda <- exp(construct_design(plotx,splineknots,p,m) %*% samps$samps)
est_lambda <- apply(sampled_lambda,1,mean)
lambda_lower <- apply(sampled_lambda,1,quantile,probs = .025)
lambda_upper <- apply(sampled_lambda,1,quantile,probs = .975)

true_lambda <- exp(truefunc(plotx))

plot(plotx,true_lambda,type='l',ylim = c(0,3),col = 'red')
lines(plotx,est_lambda)
lines(plotx,lambda_lower,lty = 'dashed')
lines(plotx,lambda_upper,lty = 'dashed')

rmse_aghq <- sqrt( mean( (est_lambda - true_lambda)^2 ) )
covr_aghq <- mean(true_lambda <= lambda_upper & true_lambda >= lambda_lower)


# Check mgcv and INLA

mgcvmod_bs <- mgcv::gam(
  y ~ s(x,bs='bs',m=c(p-1,m),k=length(splineknots)-p),
  data = dat,
  family = 'poisson'
)

mgcv_bs_pred <- predict(mgcvmod_bs,newdata = data.frame(x = plotx),se.fit = TRUE)

mgcvmod_tp <- mgcv::gam(
  y ~ s(x), # Default
  data = dat,
  family = 'poisson'
)

mgcv_tp_pred <- predict(mgcvmod_tp,newdata = data.frame(x = plotx),se.fit = TRUE)

mm = get("inla.models", inla.get.inlaEnv())
mm$latent$rw2$min.diff = NULL
assign("inla.models", mm, inla.get.inlaEnv())
inlamod <- inla(
  y ~ -1 + f(x,model = "rw2",hyper = list(theta = list(prior = "pc.prec",param = c(tmbdat$u,tmbdat$alpha)))),
  data = dat,
  family = 'poisson'
)

inlasigma <- inla.tmarginal(function(x) 1/sqrt(x),inlamod$marginals.hyperpar$`Precision for x`)
plot(inlasigma[ ,1],inlasigma[ ,2],type='l')

inlaX <- inlamod$summary.random$x$ID
inlapred <- inlamod$summary.random$x$mean
inlalower <- inlamod$summary.random$x$`0.025quant`
inlaupper <- inlamod$summary.random$x$`0.975quant`

# STAN
stanmod <- tmbstan(
  ff,
  chains = 4,
  iter = 1e04,
  warmup = 1e03,
  seed = 43780,
  algorithm = "NUTS",
  control = list(max_treedepth = 20)
)
# save(stanmod,file = paste0(globalpath,"data/stanmod-poissonbspline-20210429.RData"))


stansamps <- as.data.frame(stanmod)
# sigma
hist(exp(-stansamps$theta/2),freq=FALSE,breaks = 100)
with(logpostsigma,lines(transparam,pdf_transparam)) # AGHQ
ss <- seq(10,45,length.out = 1e03)
tt <- -2*log(ss)
with(quad,lines(ss,dnorm(tt,optresults$mode,1/sqrt(optresults$hessian))*abs(numDeriv::grad(function(x) -2*log(x),ss)),lty='dashed'))
# ks disance
numsamp <- nrow(stansamps)
quadsamp <- sample_marginal(quad,numsamp)$thetasamples[[1]]
normsamp <- rnorm(numsamp,quad$optresults$mode,1/sqrt(quad$optresults$hessian))
ks.test(stansamps$theta,quadsamp)$statistic
ks.test(stansamps$theta,normsamp)$statistic

# nice
W <- apply(stansamps[ ,grep('W',colnames(stansamps))],2,mean)
U <- W[1:ncol(P)]
beta <- W[(ncol(P)+1):length(W)]

# Construct a plot
plotx <- seq(a,b,by=0.01)
plotdat <- data.frame(x = plotx)
plotB <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,0),k=length(splineknots)-p),data = plotdat,knots = list(x = splineknots))$X
# plotX <- sparse.model.matrix(~poly(plotx,degree = m-1,raw = FALSE))
plotX <- sparse.model.matrix(~1,data = plotdat)
# Construct the plot
ploteta <- plotB %*% U + plotX %*% beta
plot(plotx,exp(ploteta),type='l',ylim=c(0,12))
# Plot the samples
samplestoplot <- stansamps[ ,sample(1:ncol(stansamps),100,replace = FALSE)]
for (i in 1:ncol(samplestoplot)) {
  WS <- samplestoplot[ ,i]
  US <- WS[1:ncol(P)]
  betaS <- WS[(ncol(P)+1):length(WS)]
  plotetaS <- plotB %*% US + plotX %*% betaS
  lines(plotx,exp(plotetaS),col = 'lightgray')
}
lines(plotx,exp(ploteta))
lines(plotx,exp(truefunc(plotx)),col='red')
points(x,y,pch=4)



# metrics
plot(plotx,exp(mgcv_bs_pred$fit),type = 'l',ylim = c(0,4))
lines(plotx,exp(mgcv_bs_pred$fit - 2*mgcv_bs_pred$se.fit),type='l',lty='dashed')
lines(plotx,exp(mgcv_bs_pred$fit + 2*mgcv_bs_pred$se.fit),type='l',lty='dashed')
lines(plotx,exp(truefunc(plotx)),col='red')

rmse_mgcv <- sqrt( mean( (exp(mgcv_bs_pred$fit) - exp(truefunc(plotx)))^2 ) )
covr_mgcv <- mean(exp(truefunc(plotx)) <= exp(mgcv_bs_pred$fit + 2*mgcv_bs_pred$se.fit) & exp(truefunc(plotx)) >= exp(mgcv_bs_pred$fit - 2*mgcv_bs_pred$se.fit))


plot(inlaX,exp(inlapred),type = 'l',ylim = c(0,4))
lines(inlaX,exp(inlalower),type='l',lty='dashed')
lines(inlaX,exp(inlaupper),type='l',lty='dashed')
lines(inlaX,exp(truefunc(inlaX)),col='red')

rmse_inla <- sqrt( mean( (exp(inlapred) - exp(truefunc(inlaX)))^2 ) )
covr_inla <- mean(exp(truefunc(inlaX)) <= exp(inlaupper) & exp(truefunc(inlaX)) >= exp(inlalower))

data.frame(
  method = c('AGHQ','MGCV','INLA'),
  rmse = c(rmse_aghq,rmse_mgcv,rmse_inla),
  coverage = c(covr_aghq,covr_mgcv,covr_inla)
)
