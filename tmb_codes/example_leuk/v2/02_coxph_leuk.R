### Coxph Leuk Example ###
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
d <- 46
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
construct_penalty <- function(x,splineknots,p,m, noise = 0.0001) {
  BD <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,m),k=length(splineknots)-p),data = data.frame(x = x),knots = list(x = splineknots))
  BB <- mgcv::smooth.construct(s(x,bs='bs',m=c(p-1,0),k=length(splineknots)-p),data = data.frame(x = x),knots = list(x = splineknots))
  # BD$S[[1]] + BB$S[[1]] # O'sullivan spline "of the third kind"
  BD$S[[1]] + diag(noise, ncol = ncol(BD$S[[1]]), nrow = nrow(BD$S[[1]])) #### add a very small noise to make the penalty matrix full rank
}


P <- as(construct_penalty(data$tpi,splineknots,p,m, noise = 0.0001),'dgTMatrix')
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
  # logPdet = sum(log(eigen(P,only.values = T)$values[1:(length(eigen(P,only.values = T)$values)-m)])),
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
compile("02_coxph_leuk.cpp")
dyn.load(dynlib("02_coxph_leuk"))

start_time <- Sys.time()
ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "02_coxph_leuk",
  silent = TRUE
)
# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)

# AGHQ
quad <- aghq::marginal_laplace_tmb(ff,7,0)

# Plot of theta posterior
logpostsigma <- compute_pdf_and_cdf(quad$marginals[[1]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)), interpolation = "spline")
with(logpostsigma,plot(transparam,pdf_transparam,type='l'))

# Inference for W
samps <- sample_marginal(quad,3000,interpolation = "spline")
end_time <- Sys.time()
runtime_aghq <- end_time - start_time
runtime_aghq

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
                  type = c(rep("Proposed",times = length(etaplotframe$x)), rep("INLA", times = length(plotINLA$x)), rep("mgcv",times = length(plotx)))  )

overall %>% filter(type != "mgcv") %>% mutate(y = exp(y)) %>% ggplot(aes(tpi,y)) + geom_line(aes(linetype = type)) + theme_classic(base_size = TEXT_SIZE)
ggsave(filename = "leukSmooth.png")




### hyper-parameter:
inla_hyper <- brinla::bri.hyperpar.plot(Inlaresult, together = F)
inla_hyper <- inla_hyper %>% filter(parameter == "SD for tpi_binned") %>% select(1:2) %>% mutate(method = "INLA")
aghq_hyper <- logpostsigma %>% select(c("transparam","pdf_transparam")) %>% mutate(method = "Proposed")
names(aghq_hyper)[1:2] <- c("x","y")
hyper <- rbind(inla_hyper,aghq_hyper)
hyper %>% ggplot(aes(x,y)) + geom_line(aes(linetype = method)) + xlim(c(0,0.5)) + xlab("SD") + ylab("density") + theme_classic(base_size = TEXT_SIZE)


theta_logprior <- function(theta, prior_alpha = tmbdat$alpha, 
                           prior_u = tmbdat$u) {
  lambda <- -log(prior_alpha)/prior_u
  log(lambda/2) - lambda * exp(-theta/2) - theta/2
}
priorfuncsigma <- function(x) (2/x) * exp(theta_logprior(-2 * log(x)))
prior <- tibble(x = hyper$x, y = priorfuncsigma(hyper$x), method = "Prior")
hyper <- rbind(prior,hyper)
hyper %>% ggplot(aes(x,y)) +  geom_line(aes(linetype = method)) + xlim(c(0,0.5)) +
  xlab("SD") + ylab("density") + 
  scale_linetype_manual(values=c("solid", "dotted","dashed")) + 
  theme_classic(base_size = TEXT_SIZE)

ggsave(filename = "leukHyper.png")








### STAN:
start_time <- Sys.time()
stanmod <- tmbstan(
  ff,
  chains = 4,
  cores = 4,
  iter = 20000,
  warmup = 10000,
  init = quad$optresults$mode,
  seed = 123
)
end_time <- Sys.time()
runtime_MCMC <- end_time - start_time
runtime_MCMC
save(stanmod, file = "stanmod.Rdat")

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



##### overall comparison:
# overall <- tibble(tpi = c(etaplotframe$x, plotINLA$x, plotx), y = c(etaplotframe$etamean,plotINLA$f,mgcv_bs_pred$fit), 
#                   type = c(rep("aghq",times = length(etaplotframe$x)), rep("INLA", times = length(plotINLA$x)), rep("mgcv",times = length(plotx))))
# overall <- rbind(overall,tibble(tpi = etaplotframe2$x, y = etaplotframe2$etamean, type = rep("MCMC",times = length(etaplotframe2$x))))
# overall %>% mutate(y = exp(y)) %>% ggplot(aes(tpi,y,color = type)) + geom_line()
# lower_frame <- etaplotframe %>% select(x,lower) %>% mutate(type = "lower")
# names(lower_frame) <- c("tpi","y", "type")
# upper_frame <- etaplotframe %>% select(x,upper) %>% mutate(type = "upper")
# names(upper_frame) <- c("tpi","y", "type")
# overall <- rbind(overall,lower_frame,upper_frame)

# overall %>% filter(type == "INLA" | type == "aghq") %>% mutate(y = exp(y)) %>% ggplot(aes(tpi,y,color = type)) + geom_line() + 
#   theme(text = element_text(size = TEXT_SIZE)) + geom_ribbon(data = etaplotframe, aes(x = x, ymin = lower, ymax = upper))

# etaplotframe %>% ggplot(aes(x = x)) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "orange", alpha = 0.4) + 
#   geom_line(aes(x = x, y = etamean, color = "Proposed")) + geom_line(data = plotINLA, aes(x = x, y = f, color = "INLA")) +
#   theme(text = element_text(size = TEXT_SIZE)) + ylab("tpi effect") + xlab("tpi") + scale_fill_manual(values = c("red","blue"))
#   
  

### KS statistics:
stansamps <- as.data.frame(stanmod)
numsamp <- nrow(stansamps)
quadsamp <- sample_marginal(quad,numsamp,interpolation = 'spline')$thetasamples[[1]]
normsamp <- rnorm(numsamp,quad$optresults$mode,1/sqrt(quad$optresults$hessian))
stansamps$sigma <- exp(-stansamps$theta/2)
ks.test(stansamps$theta,quadsamp)$statistic
ks.test(stansamps$theta,normsamp)$statistic
ks.test(stansamps$sigma,exp(-quadsamp/2))$statistic
ks.test(stansamps$sigma,exp(-normsamp/2))$statistic

# Look at the KS
# The distributions look pretty close:
hist(stansamps$theta,breaks = 50,freq=FALSE)
with(logpostsigma,lines(theta,pdf))
hist(stansamps$sigma,breaks = 50,freq=FALSE, xlim = c(0,2), ylim = c(0,5))
with(logpostsigma,lines(transparam,pdf_transparam))

ggplot(stansamps, aes(x = sigma)) + 
  geom_histogram(fill = "gray", color = "gray", stat = "density") + 
  geom_line(data = logpostsigma, aes(x = transparam, y = pdf_transparam)) + xlim(0,1.5) + theme_classic(base_size = TEXT_SIZE)
ggsave(filename = "leukHyper2.png")

# Compute the KS manually. Plot the ECDFs:
tt <- seq(-3,6,length.out=1e04)
quadecdf <- ecdf(quadsamp)(tt)
stanecdf <- ecdf(stansamps$theta)(tt)
plot(tt,quadecdf,type='l')
lines(tt,stanecdf,lty='dashed')
plotdata <- tibble(x = rep(tt,2), cdf = c(quadecdf,stanecdf), methods = rep(c("Proposed","MCMC"), each = length(tt)))
plotdata %>% ggplot(aes(x,cdf)) + geom_line(aes(linetype = methods)) + 
  theme(text = element_text(size = TEXT_SIZE))

# KS is the max absolute difference:
theKS <- max(abs(stanecdf - quadecdf))
whereistheKS <- which.max(abs(stanecdf - quadecdf))
abline(v = tt[whereistheKS])
plot(tt,abs(stanecdf - quadecdf),type='l')





