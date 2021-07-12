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
library(geostatsp, lib = lib_loc)
library(sp)
library(raster)
library(mapmisc)

TEXT_SIZE <- 25
reslist <- list(nrow = 50,ncol = 100)
### Setup_data:
data <- as.tibble(Leuk) %>% dplyr::select(c("time","cens","age","sex","wbc","tpi","xcoord","ycoord")) 
# data <- sample_n(data, 600)
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
set.seed(13579)
startingsig <- exp(-1)
startingrho <- 50*1000
# startingsig <- 1
# startingrho <- 4.22*1e04

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
compile("04_coxph_leuk.cpp")
dyn.load(dynlib("04_coxph_leuk"))


ff <- MakeADFun(data = datlist,
                parameters = paraminit,
                random = "W",
                DLL = "04_coxph_leuk",
                ADreport = FALSE,
                silent = FALSE)

# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)

timestart <- Sys.time()

quad <- aghq::marginal_laplace_tmb(
  ff,
  k = 4,
  startingvalue = c(paraminit$theta,paraminit$logkappa,paraminit$logtau)
)

save(quad,file = "quad.rda")

posterior_samples <- sample_marginal(quad,1e03)

endtime <- Sys.time()
endtime - timestart






#### RW2 Smoothing:
RWid <- n + (1:ncol(Bmat)) # Note: first index is intercept, not actually estimable here
RWsamps <- posterior_samples$samps[RWid, ]
RWpostMean <- apply(RWsamps,1,mean)
# Construct a plot
plotx <- seq(a,b,by=0.01)
plotdat <- data.frame(x = plotx)
plotB <- mgcv::smooth.construct(s(x,bs='bs',m=c(4-1,0),k=length(splineknots)-4),data = plotdat,knots = list(x = splineknots))$X
# Construct the plot
ploteta <- plotB %*% RWpostMean
samplestoplot <- RWsamps[ ,sample(1:ncol(RWsamps),200,replace = FALSE)]
plot(plotx,ploteta,type='l',ylim=c(-0.5,0.5))
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
for (i in 1:ncol(RWsamps)) {
  W <- RWsamps[ ,i]
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

etaplotframe %>% ggplot(aes(x,etamean)) + geom_line(linetype = "solid") + geom_line(aes(x,lower), linetype = "dashed") +
  geom_line(aes(x,upper), linetype = "dashed") + 
  theme_classic(base_size = TEXT_SIZE) + ylab("") + xlab("tpi")


##### Fixed effects and Hyper parameters table:
betaidx <- n + (ncol(Bmat)+1):(ncol(BX)) # Note: first index is intercept, not actually estimable here
betasamps <- posterior_samples$samps[betaidx, ]
betamean <- apply(betasamps,1,mean)
betasd <- apply(betasamps,1,sd)
beta2.5 <- apply(betasamps,1,quantile,probs = .025)
beta97.5 <- apply(betasamps,1,quantile,probs = .975)





########## Spatial effect:
# Simulate the spatial fields
simulate_spatial_fields <- function(U,
                                    theta,
                                    pointsdata,
                                    resolution = list(nrow = 100,ncol = 100)) {
  # U: matrix of samples, each column is a sample
  # theta: tibble of theta values
  # Draw from U*|U
  fieldlist <- vector(mode = 'list',length = nrow(theta))
  for (i in 1:length(fieldlist)) {
    fielddat <- pointsdata
    fielddat@data <- data.frame(w = as.numeric(U[ ,i]))
    # Back-transform the Matern params
    kappa <- exp(theta$theta2[i])
    tau <- exp(theta$theta3[i])
    sig <- get_sigma(kappa, tau)
    rho <- get_rho(kappa, tau)
    # Simulate from the two fields
    capture.output({
      fieldlist[[i]] <- geostatsp::RFsimulate(
        model = c("variance" = sig^2,"range" = rho,"shape" = maternconstants$nu),
        data = fielddat,
        x = raster(fielddat,nrows = resolution$nrow,ncols = resolution$ncol),
        n = 1
      )
    })
  }
  brick(fieldlist)
}
simstodo <- sample(ncol(posterior_samples$samps),100,replace = FALSE)
fieldsbrick <- simulate_spatial_fields(
  U = posterior_samples$samps[1:n,simstodo],
  theta = posterior_samples$theta[simstodo, 2:3],
  pointsdata = pointsdata,
  resolution = list(nrow = 400,ncol = 200)
)
save(fieldsbrick, file = "fieldsbrick.rda")


###### Spatial plot:
ukBorderLL = raster::getData("GADM", country='GBR', level=3) # Regions
ukBorder = spTransform(ukBorderLL, projection(pointsdata))
# ukBorder = ukBorder[ukBorder$NAME_1 %in% c("England","Wales"), ]
# ukBorder = raster::crop(ukBorder, extent(pointsdata))
# TODO: Plot only polygons that have a point in them, this isn't quite what's being done
pointsinpoly <- pointsdata %over% ukBorder
pointsinpolyID <- unique(pointsinpoly$GID_2)
ukBorder <- ukBorder[ukBorder$GID_2 %in% pointsinpolyID, ]
# Get the outer border
ukBorderouter <- rgeos::gUnaryUnion(ukBorder)

simfieldsmean <- mean(exp(fieldsbrick))
simfieldsexceedence <- mean(fieldsbrick > log(1.5))

# MEAN
plotraster <- simfieldsmean

predcols <- mapmisc::colourScale(
  plotraster,
  breaks = quantile(values(plotraster),probs = (0:9)/9),
  style = "fixed",
  col = "Spectral",
  rev = TRUE,
  # transform='log',
  dec = -log10(0.05)
)

colvals <- 100
bigvalues <- quantile(values(plotraster),probs = (0:(colvals-1))/(colvals-1))
plotraster <- mask(plotraster, ukBorder)

mapmisc::map.new(pointsdata)
plot(plotraster,
     col = predcols$col,
     breaks = predcols$breaks,
     legend=FALSE, add=TRUE)
plot(plotraster, col=predcols$col, breaks=predcols$breaks, legend=FALSE, add=TRUE)
plot(ukBorder, add=TRUE,border=mapmisc::col2html("black", 0.5), lwd=0.5)
plot(ukBorderouter,add = TRUE)
points(pointsdata,pch = ".")
mapmisc::legendBreaks('topright', predcols, cex=1.5, bty='n',inset=0)



# EXCEEDENCE PROBABILITIES
plotraster <- simfieldsexceedence

predcols <- mapmisc::colourScale(
  plotraster,
  breaks = c(0,0.05,0.1,0.15,0.2,0.3,0.5,0.65),
  style = "fixed",
  col = "Spectral",
  rev = TRUE
)

plotraster <- mask(plotraster, ukBorder)

mapmisc::map.new(pointsdata)
plot(plotraster,
     col = predcols$col,
     breaks = predcols$breaks,
     legend=FALSE, add=TRUE)
plot(ukBorder, add=TRUE,border=mapmisc::col2html("black", 0.5), lwd=0.5)
plot(ukBorderouter,add = TRUE)
points(pointsdata,pch = ".")
mapmisc::legendBreaks('topright', predcols, cex=1.5, bty='n')




































#### Study the factorization time of H matrix:
Mode <- quad$modesandhessians$mode[[1]]
Theta_ori <- c(quad$modesandhessians$theta2[1], quad$modesandhessians$theta3[1], quad$modesandhessians$theta1[1])
Theta <- Theta_ori
Theta[1] <- log(get_sigma(exp(Theta_ori[1]), exp(Theta_ori[2])))
Theta[2] <- log(get_rho(exp(Theta_ori[1]), exp(Theta_ori[2])))
Theta[3] <- -0.5 * Theta_ori[3]


##### From the novel method (without epsilon noise):
Hnew <- quad$modesandhessians$H[[1]]
BeginTime <- Sys.time()
chol(Hnew)
EndTime <- Sys.time()
newTime <- EndTime - BeginTime


##### From the old method (with epsilon noise):

### Compute the old C matrix:

# Setup the fixed effect term and the design matrix
D <- abcoxph:::create_diff_matrix(n)
X <- as(sparse.model.matrix(times ~ -1 + age + sex + wbc ,data = data),'dgTMatrix')
BX <- as(cbind(B,X),'dgTMatrix')

Amat <- Diagonal(n = nrow(data),x = 1)
censor <- data$cens[-1] # 1 == not censored, confusing but more useful.

# Zmat is the differenced design matrix
Zmat <- D %*% cbind(BX,Amat) 

make_delta <- function(W) {
  as.numeric(Zmat %*% cbind(W))
}

Delta <- make_delta(Mode)
oldW <- c(Delta, Mode)

######### some necessary functions
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

oldC <- hessian_log_likelihood(Mode)
oldC <- bdiag(oldC, diag(0,length(Mode),length(Mode)))

### Compute the Old Q matrix:

Q_matrix <- function(theta) {
  tau <- exp(7)
  I <- Diagonal(n-1)
  # theta = log(sigma), log(rho), log(Sigma)
  theta <- as.numeric(unname(theta))
  # RW2 is P matrix
  PrecRW <-  exp(theta[3]) * P
  # Matern
  PrecMatern <- geostatsp::matern(
    pointsdata,
    param = c("variance" = exp(2 * theta[1]),"range" = exp(theta[2]),"shape" = 1),
    type = "precision"
  )
  # fixed effect
  PrecBeta <- beta_prec * diag(p)
  # Transformed Design matrix:
  DA <- D %*% Amat
  DB <- D %*% Bmat
  DX <- D %*% Xmat
  tau * rbind(
    cbind(I, -DA, -DB, -DX),
    cbind(-t(DA), (1/tau)* PrecMatern + t(DA)%*%DA, t(DA)%*%DB, t(DA)%*%DX),
    cbind(-t(DB), t(DB)%*%DA, (1/tau)*PrecRW + t(DB)%*%DB, t(DB)%*%DX),
    cbind(-t(DX), t(DX)%*%DA, t(DX)%*%DB, (1/tau)*PrecBeta + t(DX)%*%DX)
  )
}

oldQ <- Q_matrix(Theta)

### Factorize an old H matrix:
oldH <- oldC + oldQ
BeginTime <- Sys.time()
chol(oldH)
EndTime <- Sys.time()
oldTime <- EndTime - BeginTime

#### Compare factorization time between two methods:
newTime
oldTime

#### Compare the memory difference between two methods:
print(object.size(Hnew), units = "Mb")
print(object.size(oldH), units = "Mb")
