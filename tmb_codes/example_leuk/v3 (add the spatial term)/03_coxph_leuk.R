### Coxph Leuk Example ###
library(tidyverse)
library(aghq)
library(mgcv)
library(Matrix)
library(TMB)
library(INLA)
library(tmbstan)
library(brinla)
library(geostatsp)
library(sp)
library(raster)
library(mapmisc)

precompile()
TEXT_SIZE = 25
### sigma for Spatial effect
sigma_u <- 1
sigma_alpha <- .5
rho_u <- 20 * 1000
rho_alpha <- .5
### Tau_RW for RW2 effect
Tau_RW_u <- 2
Tau_RW_alpha <- 0.5
### For the fixed effect
beta_prec <- .001

### Setup_data:
set.seed(1234)
data <- Leuk %>% select(c("time","cens","age","sex","wbc","tpi","xcoord","ycoord"))
names(data)[1] <- "times"
data <- abcoxph:::arrange_data(data)
data$ranks <- rank(data$times, ties.method = "min")
PROJTOUSE <- mapmisc::omerc(c(-3.055,53.365),angle=0)
crstouse <- CRS("+init=epsg:27700")
pointsdata <- SpatialPointsDataFrame(
  coords = 90*1000*dplyr::select(data,xcoord,ycoord),
  data = dplyr::select(data,-xcoord,-ycoord),
  proj4string = PROJTOUSE
)
pointsdata <- spTransform(pointsdata,crstouse)


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

# Setup the RW2 term:
P <- as(construct_penalty(data$tpi,splineknots,p,m, noise = 0.0001),'dgTMatrix')
B <- as(construct_design(data$tpi,splineknots,p,m),'dgTMatrix')




# Setup the fixed effect term and the design matrix
D <- abcoxph:::create_diff_matrix(n)
X <- as(sparse.model.matrix(times ~ -1 + age + sex + wbc ,data = data),'dgTMatrix')
BX <- as(cbind(B,X),'dgTMatrix')

Amat <- Diagonal(n = nrow(data),x = 1)
censor <- data$cens[-1] # 1 == not censored, confusing but more useful.

# Zmat is the differenced design matrix
Zmat <- D %*% BX

make_delta <- function(W) {
  as.numeric(Zmat %*% cbind(W))
}

# Dimensions
p <- ncol(Xmat)
d <- ncol(Amat)
Wd <- ncol(Zmat)
stopifnot(p+d == Wd)


# Likelihood and derivatives

compute_one_denominator <- function(delta,i) {
  # All of the likelihood quantities require that denominator
  # vector for each observation. It's a cumulative sum. Write
  # one function that computes it efficiently.
  # delta: vector of length n
  # i: index of denominator you want
  n <- length(delta)
  dd <- delta[i] - delta[i:n]
  exp(matrixStats::logSumExp(dd)) - 1
}

compute_denominator <- function(delta) {
  map(1:length(delta),~compute_one_denominator(delta,.x)) %>% reduce(c)
}

log_likelihood <- function(W) {
  delta <- make_delta(W)
  denom <- compute_denominator(delta)
  -sum(censor * log(1 + denom))
}

grad_log_likelihood_one <- function(W,i) {
  delta <- make_delta(W)
  n <- length(delta)
  if (censor[i] == 0) return(sparseVector(0,0,n))
  denom <- compute_one_denominator(delta,i)
  dd <- delta[i] - delta[i:n]
  out <- exp(dd) / (1 + denom)
  out <- c(rep(0,i-1),out)
  out[i] <- out[i] - 1
  out
}

grad_log_likelihood_subset <- function(W,I) {
  # I: subset of 1...n
  map(I,~grad_log_likelihood_one(W,.x)) %>% reduce(~.x + .y)
}

grad_log_likelihood <- function(W) grad_log_likelihood_subset(W,which(censor==1))


make_hess_vec <- function(delta,i) {
  # Make the vector that is used to create the hessian
  n <- length(delta)
  denom <- compute_one_denominator(delta,i)
  if (censor[i] == 0) return(sparseVector(0,0,n))
  dd <- delta[i] - delta[i:n]
  out <- exp(dd) / (1 + denom)
  c(rep(0,i-1),out)
}


hessian_log_likelihood <- function(W) {
  delta <- make_delta(W)
  n <- length(delta)
  gg <- map(which(censor==1),~make_hess_vec(delta,.x))
  diag(as.numeric(reduce(gg,~.x+.y))) - tcrossprod(gg %>% reduce(cbind))
}





# Prior
Q_matrix <- function(theta) {
  # theta = log(sigma), log(rho), log(Sigma)
  theta <- as.numeric(unname(theta))
  # RW2 is P matrix
  P2 <-  exp(theta[3]) * P
  # Matern
  mm <- geostatsp::matern(
    pointsdata,
    param = c("variance" = exp(2 * theta[1]),"range" = exp(theta[2]),"shape" = 1),
    type = "precision"
  )
  # fixed effect
  bb <- beta_prec * diag(p)
  rbind(
    cbind(P2,Matrix(0,nrow = nrow(P2),ncol = (ncol(mm)+ncol(X)),sparse = FALSE)),
    cbind(Matrix(0,nrow = p,ncol = ncol(P2),sparse = FALSE),bb, Matrix(0,nrow = p,ncol = ncol(mm),sparse = FALSE)),
    cbind(Matrix(0,nrow = nrow(mm),ncol = (ncol(P2) + ncol(X)),sparse = FALSE),mm)
  )
}


logsigmalogprior <- function(theta,sigma0,prior_alpha) {
  # theta = log(sigma)
  alpha2 <- prior_alpha
  lambda2 <- -log(alpha2) / sigma0
  log(lambda2) - exp(theta) * lambda2 + theta
}

logTaulogprior <- function(theta,Tau_RW_u,prior_alpha) {
  # theta is now log precision
  alpha2 <- prior_alpha
  lambda2 <- -log(alpha2) / Tau_RW_u
  log(lambda2) - exp(-theta/2) * lambda2 - theta/2
}

logrhologprior <- function(theta,rho0,prior_alpha) {
  # theta = log(rho)
  d <- 2
  alpha1 <- prior_alpha
  lambda1 <- -log(alpha1) * rho0^(d/2)
  log(d/2) + log(lambda1) - (d/2)*theta - lambda1 * exp(-theta * d/2)
}

logprior_theta <- function(theta) {
    logsigmalogprior(theta[1],sigma0 = sigma_u,prior_alpha = sigma_alpha) + 
    logrhologprior(theta[2],rho0 = rho_u,prior_alpha = rho_alpha) +
    logTaulogprior(theta[3],Tau_RW_u = Tau_RW_u,prior_alpha = Tau_RW_alpha)
}

log_prior_W <- function(W,theta,Q = NULL) {
  if (is.null(Q)) Q <- Q_matrix(theta)
  -(1/2) * as.numeric(crossprod(W,crossprod(Q,W))) +(1/2)*as.numeric(determinant(Q,logarithm = TRUE)$modulus)
}