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


reslist <- list(nrow = 50,ncol = 100)
### Setup_data:
data <- as.tibble(Leuk) %>% dplyr::select(c("time","cens","age","sex","wbc","tpi","xcoord","ycoord"))
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



# Design matrices
Amat <- Diagonal(nrow(data))
Bmat <- as(construct_design(data$tpi,splineknots,p,m),'dgTMatrix')
Xmat <- as.matrix(data[,3:5])
BX = cbind(Bmat,Xmat)
# Design matrix: zip model and risk model are the same
design <- bdiag(
  # ZIP
  cbind(
    Amat,
    BX
  )
)



## Dimensions
n <- nrow(Xmat) # Number of obs
p <- ncol(Xmat) # Number of betas
m <- ncol(Amat) # Number of spatial points
Wd <- ncol(design) # Number of total params
D <- abcoxph:::create_diff_matrix(n)

## Prior distributions ##
# Use the same prior for both sets of Matern params
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

# PC Prior for kappa,tau
maternconstants <- list()
maternconstants$d <- 2 # Dimension of spatial field, fixed
maternconstants$nu <- 1 # Shape param, fixed
get_kappa <- function(sigma,rho) sqrt(8*maternconstants$nu)/rho
get_tau <- function(sigma,rho) sigma * get_kappa(sigma,rho)^(maternconstants$nu) * sqrt(gamma(maternconstants$nu + maternconstants$d/2) * (4*pi)^(maternconstants$d/2) / gamma(maternconstants$nu))
get_sigma <- function(kappa,tau) tau / (kappa^(maternconstants$nu) * sqrt(gamma(maternconstants$nu + maternconstants$d/2) * (4*pi)^(maternconstants$d/2) / gamma(maternconstants$nu)))
get_rho <- function(kappa,tau) sqrt(8*maternconstants$nu) / kappa




datlist <- list(
  # Response
  ranks = as.integer(data$ranks),
  cens = as.integer(data$cens),
  design = design,
  # Penalty matrix
  P = P,
  # Differencing matrix
  D = D,
  # Log determinant of penalty matrix (without the sigma part)
  logPdet = as.numeric(determinant(P,logarithm = TRUE)$modulus),
  nu = maternconstants$nu,
  rho_u = rho_u,
  rho_alpha = rho_alpha,
  sigma_u = sigma_u,
  sigma_alpha = sigma_alpha,
  # Prior params for RW2 term
  u = 2,
  alpha = 0.5,
  DS = raster::pointDistance(pointsdata,lonlat = FALSE),
  betaprec = beta_prec
)

# NOTE: for some initial values of W, TMB's inner optimization seems to fail
# This was tried over a bunch of random initialization and most worked, and all
# gave the same optimum. But this is why we set the seed here and use a random start.
set.seed(4564)
# startingsig <- exp(-1)
# startingrho <- 50*1000
startingsig <- 1
startingrho <- 4.22*1e04

paraminit <- list(
  W = rnorm(ncol(design)),
  theta = 0, # for the RW2 term
  # W = rep(0, ncol(design)),
  logkappa = log(get_kappa(startingsig,startingrho)),
  logtau = log(get_tau(startingsig,startingrho))
  # logkappa = -10,
  # logtau = -2
)


# Compile TMB template-- only need to do once
# compile("04_coxph_leuk.cpp")
dyn.load(dynlib("04_coxph_leuk"))


ff <- MakeADFun(data = datlist,
                parameters = paraminit,
                random = "W",
                DLL = "04_coxph_leuk",
                ADreport = TRUE,
                silent = TRUE)

quad <- aghq::marginal_laplace_tmb(
  ff,
  1,
  startingvalue = c(0,0,0)
)


##### current problem: the initial values seem to be not finite...









