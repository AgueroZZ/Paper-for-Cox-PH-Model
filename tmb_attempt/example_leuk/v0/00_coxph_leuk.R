### Coxph Leuk Example ###
library(tidyverse)
library(aghq)
library(mgcv)
library(Matrix)
library(TMB)
library(INLA)
library(tmbstan)
library(brinla)
library(INLA)

precompile()
TEXT_SIZE = 25

### Setup_data:
set.seed(1234)
data <- Leuk %>% select(c("time","cens","age","sex","wbc","tpi"))
names(data)[1] <- "times"
data <- abcoxph:::arrange_data(data)
data$ranks <- rank(data$times, ties.method = "min")
hist(data$times,breaks = 100)

## setup smoothing part:
a <- min(data$tpi)
b <- max(data$tpi) # boundary
n <- nrow(data)
# Order of spline
p <- 4 # 4 = cubic
# Order of derivative penalty
m <- 2
# Number of INTERIOR knots
d <- 42
# Number of knots
T <- d + p
# The knots
intknots <- seq(a,b,length.out = d)
leftknots <- seq(min(intknots)-(p-1),min(intknots)-1,by=1)
rightknots <- seq(max(intknots)+1,max(intknots)+p-1,by=1)
splineknots <- sort(unique(c(leftknots,intknots,rightknots)))
construct_design <- function(x,splineknots,p,m) {
  BB <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,0),k=length(splineknots)-p),data = data.frame(x = x),knots = list(x = splineknots))
  B <- BB$X
  B <- as(B,"dgTMatrix")
  B
}
construct_penalty <- function(x,splineknots,p,m) {
  BD <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,m),k=length(splineknots)-p),data = data.frame(x = x),knots = list(x = splineknots))
  BB <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,0),k=length(splineknots)-p),data = data.frame(x = x),knots = list(x = splineknots))
  BD$S[[1]] + BB$S[[1]] # O'sullivan spline "of the third kind"
}


P <- as(construct_penalty(data$tpi,splineknots,p,m),'dgTMatrix')
B <- as(construct_design(data$tpi,splineknots,p,m),'dgTMatrix')
D <- abcoxph:::create_diff_matrix(n)
X <- as.matrix(data[,3:5])
BX <- as(cbind(B,X),'dgTMatrix')

### Setup TMB:
tmbdat <- list(
  # Design matrix
  BX = BX,
  # Penalty matrix
  P = P,
  # Differencing matrix
  D = D,
  # Log determinant of penalty matrix (without the sigma part)
  logPdet = as.numeric(determinant(P,logarithm = TRUE)$modulus),
  # Response
  ranks = as.integer(data$ranks),
  cens = as.integer(data$cens),
  # Prior params
  u = 2,
  alpha = 0.5
)

tmbparams <- list(
  W = rep(0,ncol(BX)), # W = c(U); U = B-Spline coefficients
  theta = 0 # -2log(sigma)
)

# TMB function template
compile("00_coxph_leuk.cpp")
dyn.load(dynlib("00_coxph_leuk"))

ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "00_coxph_leuk",
  silent = TRUE
)
# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)

# AGHQ
quad <- aghq::marginal_laplace_tmb(ff,15,0)

# Plot of theta posterior
logpostsigma <- compute_pdf_and_cdf(quad$marginals[[1]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)), interpolation = "spline")
with(logpostsigma,plot(transparam,pdf_transparam,type='l'))

# Inference for W
samps <- sample_marginal(quad,3000,interpolation = "spline")

# Plot of curve
# Posterior mean
W <- apply(samps$samps,1,mean)
U <- W[1:ncol(P)]
beta <- W[(ncol(P)+1):ncol(BX)]


# Construct a plot
plotx <- seq(a,b,by=0.01)
plotdat <- data.frame(x = plotx)
plotB <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,0),k=length(splineknots)-p),data = plotdat,knots = list(x = splineknots))$X
# Construct the plot
ploteta <- plotB %*% U
ploteta <- ploteta - mean(ploteta)
plot(plotx,ploteta,type='l',ylim=c(-0.5,0.5))
# Plot the samples
samplestoplot <- samps$samps[ ,sample(1:ncol(samps$samps),200,replace = FALSE)]
for (i in 1:ncol(samplestoplot)) {
  WS <- samplestoplot[ ,i]
  US <- WS[1:ncol(P)]
  plotetaS <- plotB %*% US
  plotetaS <- plotetaS - mean(plotetaS)
  lines(plotx,plotetaS,col = 'lightgray')
}
lines(plotx,ploteta)


# get pointwise SD of eta
etaplot <- list()
for (i in 1:ncol(samps$samps)) {
  W <- samps$samps[ ,i]
  U <- W[1:ncol(P)]
  eta <- as.numeric(plotB %*% U)
  etaplot[[i]] <- data.frame(
    x = plotx,
    eta = eta - mean(eta)
  )
}
etaplotframe <- purrr::reduce(etaplot,rbind) %>%
  group_by(x) %>%
  summarize(etamean = mean(eta),etasd = sd(eta)) %>%
  mutate(lower = etamean - 2*etasd,upper = etamean + 2*etasd)
with(etaplotframe,plot(x,etamean,type='l',ylim = c(-0.5,0.5)))
with(etaplotframe,lines(x,lower,type='l',lty='dashed'))
with(etaplotframe,lines(x,upper,type='l',lty='dashed'))





# Check mgcv
mgcv_knots <- data.frame(x = splineknots)
mgcvmod_bs <- mgcv::gam(
  times ~ age + sex + wbc + s(tpi,bs='bs',m=c(p-1,m),k=length(splineknots)-p),
  data = data,
  family = cox.ph(),
  knots = mgcv_knots,
  weights = cens
)
mgcv_bs_pred <- predict(mgcvmod_bs,newdata = data.frame(tpi = plotx, age = rep(0, length(plotx)), sex = rep(0, length(plotx)), wbc = rep(0, length(plotx))),se.fit = TRUE)
mgcv_bs_pred$fit <- mgcv_bs_pred$fit - mean(mgcv_bs_pred$fit)
plot(plotx,mgcv_bs_pred$fit,type = 'l',ylim = c(-0.5,0.5))
lines(plotx,(mgcv_bs_pred$fit - 2*mgcv_bs_pred$se.fit),type='l',lty='dashed')
lines(plotx,(mgcv_bs_pred$fit + 2*mgcv_bs_pred$se.fit),type='l',lty='dashed')







### INLA:
prior.prec <- list(prec = list(prior = "pc.prec",
                               param = c(tmbdat$u, tmbdat$a)))
data$tpi_binned <- abcoxph:::bin_covariate(data$tpi,bins = 50, type = "equal")
formula <- inla.surv(times,cens) ~ age + sex + wbc + f(tpi_binned,model = 'rw2',constr = T, hyper = prior.prec)
Inlaresult <- inla(formula = formula,data = data, family = "coxph")
fhat <- Inlaresult$summary.random$tpi_binned$mean
fup <- Inlaresult$summary.random$tpi_binned$`0.975quant`
flo <- Inlaresult$summary.random$tpi_binned$`0.025quant`
plotINLA <- data.frame(x = Inlaresult$summary.random$tpi_binned$ID, f = fhat, up = fup, lo = flo)
plot(plotINLA$x,plotINLA$f,type = 'l',ylim = c(-0.5,0.5))
lines(plotINLA$x,plotINLA$up,type='l',lty='dashed')
lines(plotINLA$x,plotINLA$lo,type='l',lty='dashed')
### INLA's inferred baseline
INLA_base <- bri.basehaz.plot(Inlaresult, plot = F)
INLA_base$basehaz <- exp(INLA_base$basehaz)
INLA_base <- data.frame(INLA_base)
ggplot(INLA_base, aes(time,basehaz)) + geom_line()







##### overall comparison:
overall <- tibble(tpi = c(etaplotframe$x, plotINLA$x, plotx), y = c(etaplotframe$etamean,plotINLA$f,mgcv_bs_pred$fit), 
                  type = c(rep("aghq",times = length(etaplotframe$x)), rep("INLA", times = length(plotINLA$x)), rep("mgcv",times = length(plotx)))  )

overall %>% mutate(y = exp(y)) %>% ggplot(aes(tpi,y,color = type)) + geom_line()






### STAN:
start_time <- Sys.time()
stanmod <- tmbstan(
  ff,
  chains = 4,
  cores = 4,
  iter = 40000,
  warmup = 30000,
  init = quad$optresults$mode,
  seed = 123
)
end_time <- Sys.time()
runtime_MCMC <- end_time - start_time

summ <- summary(stanmod)
U_mcmc <- summ$summary[1:(d+m),1]
ploteta_mcmc <- plotB %*% U_mcmc
ploteta_mcmc <- ploteta_mcmc - mean(ploteta_mcmc)
plot(plotx,ploteta_mcmc,type='l',ylim=c(-0.5,0.5))


# Construct the plot
ploteta_mcmc <- plotB %*% U_mcmc
ploteta_mcmc <- ploteta_mcmc - mean(ploteta_mcmc)
plot(plotx,ploteta_mcmc,type='l',ylim=c(-0.5,0.5))
# Plot the samples
samps_mcmc <- extract(stanmod)
samps_mcmc <- samps_mcmc$W
samplestoplot <- samps_mcmc[sample(1:nrow(samps_mcmc),1000,replace = FALSE),]
etaplot2 <- list()
for (i in 1:ncol(samplestoplot)) {
  WS <- samplestoplot[i,]
  US <- WS[1:ncol(P)]
  plotetaS <- plotB %*% US
  plotetaS <- plotetaS - mean(plotetaS)
  lines(plotx,plotetaS,col = 'lightgray')
  etaplot2[[i]] <- data.frame(
    x = plotx,
    eta = plotetaS
  )
}
lines(plotx,ploteta)

etaplotframe2 <- purrr::reduce(etaplot2,rbind) %>%
  group_by(x) %>%
  summarize(etamean = mean(eta),etasd = sd(eta)) %>%
  mutate(lower = etamean - 2*etasd,upper = etamean + 2*etasd)
with(etaplotframe2,plot(x,etamean,type='l',ylim = c(-0.5,0.5)))
with(etaplotframe2,lines(x,lower,type='l',lty='dashed'))
with(etaplotframe2,lines(x,upper,type='l',lty='dashed'))




