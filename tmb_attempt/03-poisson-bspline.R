### Poisson B-Spline regression ###

library(tidyverse)
library(aghq)
library(mgcv)
library(Matrix)
library(TMB)
precompile()
library(INLA)


set.seed(8079)
n <- 1000 # sample size
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
BD <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,m),k=length(splineknots)-p),data = dat,knots = list(x = splineknots))
BB <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,0),k=length(splineknots)-p),data = dat,knots = list(x = splineknots))

# rankMatrix(BD$S[[1]]) # 10
# rankMatrix(BB$S[[1]]) # 12
# rankMatrix(BD$S[[1]] + BB$S[[1]]) # 12

P <- BD$S[[1]] + BB$S[[1]] # O'sullivan spline "of the third kind"
P <- as(P,"dgTMatrix")
B <- BB$X
B <- as(B,"dgTMatrix")


# Polynomial basis
# degree m-1 (order m) orthognal polynomial
# X <- sparse.model.matrix(~poly(x,degree = m-1,raw = FALSE))
X <- sparse.model.matrix(~1,data = dat)



tmbdat <- list(
  # Design matrix
  BX = as(cbind(B,X),'dgTMatrix'),
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
  W = rep(0,ncol(B) + ncol(X)), # W = c(U,beta); U = B-Spline coefficients, beta = polynomial coefficients starting with intercept
  theta = 0 # -2log(sigma)
)

# TMB function template
compile("03_poisson_bspline.cpp")
dyn.load(dynlib("03_poisson_bspline"))

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

# get pointwise SD of eta
etaplot <- list()
for (i in 1:ncol(samps$samps)) {
  W <- samps$samps[ ,i]
  U <- W[1:ncol(P)]
  beta <- W[(ncol(P)+1):length(W)]
  eta <- as.numeric(plotB %*% U + plotX %*% beta)
  etaplot[[i]] <- data.frame(
    x = plotx,
    eta = eta
  )
}
etaplotframe <- purrr::reduce(etaplot,rbind) %>%
  group_by(x) %>%
  summarize(etamean = mean(eta),etasd = sd(eta)) %>%
  mutate(lower = etamean - 2*etasd,upper = etamean + 2*etasd)

with(etaplotframe,plot(x,etamean,type='l',ylim = c(-10,10)))
with(etaplotframe,lines(x,lower,type='l',lty='dashed'))
with(etaplotframe,lines(x,upper,type='l',lty='dashed'))


# Note: smoothing gets better with larger n but smoothing param posterior seems to change range with n. Hm.
# Also, the samples are crazy, so it's still wrong. But it runs.
## Next: add the rest of the quad points. Then update aghq!

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

inlapred <- inlamod$summary.random


# plot
plot(plotx,exp(mgcv_bs_pred$fit),type = 'l',ylim = c(0,4))
lines(plotx,exp(mgcv_bs_pred$fit - 2*mgcv_bs_pred$se.fit),type='l',lty='dashed')
lines(plotx,exp(mgcv_bs_pred$fit + 2*mgcv_bs_pred$se.fit),type='l',lty='dashed')
lines(plotx,exp(truefunc(plotx)),col='red')

