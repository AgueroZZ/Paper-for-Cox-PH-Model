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
Zmat <- D %*% cbind(BX,Amat) 

make_delta <- function(W) {
  as.numeric(Zmat %*% cbind(W))
}

# Dimensions
p <- ncol(X)
d <- ncol(Amat)
Wd <- ncol(Zmat)


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


#### Unnormalized posterior 

log_posterior_W <- function(W,theta,Q = NULL) {
  if (is.null(Q)) Q <- Q_matrix(theta)
  log_prior_W(W,theta,Q) + log_likelihood(W) + logprior_theta(theta)
}

grad_log_posterior_W <- function(W,theta,Q = NULL) {
  if (is.null(Q)) Q <- Q_matrix(theta)
  as.numeric(-Q %*% cbind(W) + t(Zmat) %*% grad_log_likelihood(W))
}

H_matrix <- function(W,theta) {
  C <- hessian_log_likelihood(W)
  Q <- Q_matrix(theta)
  Q + crossprod(Zmat,crossprod(C,Zmat))
}




## Fit model ----
ff <- list(
  fn = log_posterior_W,
  gr = grad_log_posterior_W,
  he = function(W,theta) -1 * H_matrix(W,theta)
)
tm <- Sys.time()
cat("Fitting model, time = ",format(tm),"\n")

coxphmod <- marginal_laplace(
  ff,
  k = 3,
  # theta starting values from Brown (2015), geostatsp/diseasemapping software paper
  startingvalue = list(W = rep(0,Wd),theta = c(-1,log(50*1000),0)),
  control = list(method = 'BFGS',inner_method = 'trust', negate = F)
)

cat("Finished model, took ",format(difftime(Sys.time(),tm,units = 'secs')),"\n")




################ Posterior summary:
set.seed(123)
posterior_samples <- sample_marginal(coxphmod,1e03)


# Plot of hyper-parameters:
sigmaprior <- function(sigma,u,alpha) dexp(sigma,-log(alpha)/u) # Same for tau as well
rhoprior <- function(rho,u,alpha) (1/rho^2) * dexp(1/rho,-log(alpha) * u)
sigmapdf <- compute_pdf_and_cdf(coxphmod$marginals[[1]],list(totheta = log,fromtheta = exp),seq(-10,-.2,by=.01))
rhopdf <- compute_pdf_and_cdf(coxphmod$marginals[[2]],list(totheta = log,fromtheta = exp),seq(0,12.5,by=.01))
Sigmapdf <- compute_pdf_and_cdf(coxphmod$marginals[[3]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)), interpolation = "spline")


### Hyper 1:
sigmaplot <- sigmapdf %>%
  ggplot(aes(x = transparam,y = pdf_transparam)) +
  theme_classic() +
  geom_line() +
  stat_function(fun = sigmaprior,args = list(u = sigma_u,alpha = sigma_alpha),linetype = 'dashed') +
  # labs(x = expression(sigma),y = "Density") +
  labs(x = "",y = "") +
  theme(text = element_text(size = TEXT_SIZE))
### Hyper 2:
rhoplot <- rhopdf %>%
  ggplot(aes(x = transparam,y = 1e05*pdf_transparam)) +
  theme_classic() +
  geom_line() +
  geom_line(data = tibble(x = seq(0.1,2e05,length.out = 1e05),y = 1e05*rhoprior(x,rho_u,rho_alpha)),aes(x = x,y = y),linetype = "dashed") +
  scale_x_continuous(breaks = seq(0,1e06,by = 2.5e04),labels = function(x) x/1000) +
  theme(text = element_text(size = TEXT_SIZE)) +
  coord_cartesian(xlim = c(0,1.5e05)) +
  labs(x = "",y = "")
### Hyper 3:
Sigmaplot <- Sigmapdf %>%
  ggplot(aes(x = transparam,y = pdf_transparam)) +
  theme_classic() +
  geom_line() +
  stat_function(fun = sigmaprior,args = list(u = Tau_RW_u,alpha = Tau_RW_alpha),linetype = 'dashed') +
  # labs(x = expression(sigma),y = "Density") +
  labs(x = "",y = "") +
  theme(text = element_text(size = TEXT_SIZE))



#### RW2 Smoothing:
RWid <- 1:ncol(B) # Note: first index is intercept, not actually estimable here
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



##### Fixed effects and Hyper parameters table:
betaidx <- (ncol(B)+1):(ncol(BX)) # Note: first index is intercept, not actually estimable here
betasamps <- posterior_samples$samps[betaidx, ]

betamean <- apply(betasamps,1,mean)
betasd <- apply(betasamps,1,sd)
beta2.5 <- apply(betasamps,1,quantile,probs = .025)
beta97.5 <- apply(betasamps,1,quantile,probs = .975)

# Sigmas
sigma_and_rho_means <- compute_moment(coxphmod$normalized_posterior,exp)
sigma_and_rho_sd <- sqrt(compute_moment(coxphmod$normalized_posterior,function(x) ((exp(x) - sigma_and_rho_means)^2)))
sigma_quants <- exp(compute_quantiles(coxphmod$marginals[[1]]))
rho_quants <- exp(compute_quantiles(coxphmod$marginals[[2]]))
Sigma_quants <- exp(-0.5*compute_quantiles(coxphmod$marginals[[3]]))
makelogPrectoSig <- function(x) {exp(-0.5*x)}
sigma_and_rho_means[3] <- compute_moment(coxphmod$normalized_posterior, makelogPrectoSig)[3]
sigma_and_rho_sd[3] <- sqrt(compute_moment(coxphmod$normalized_posterior,function(x) ((exp(-0.5*x) - sigma_and_rho_means)^2)))[3]

coeftable <- tibble(
  mean = c(betamean,sigma_and_rho_means),
  sd = c(betasd,sigma_and_rho_sd),
  q2.5 = c(beta2.5,sigma_quants[1],rho_quants[1],Sigma_quants[1]),
  q97.5 = c(beta97.5,sigma_quants[2],rho_quants[2],Sigma_quants[2])
)



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
    sig <- exp(theta$theta1[i])
    rho <- exp(theta$theta2[i])
    # Simulate from the two fields
    capture.output({
      fieldlist[[i]] <- geostatsp::RFsimulate(
        model = c("variance" = sig^2,"range" = rho,"shape" = 2),
        data = fielddat,
        x = raster(fielddat,nrows = resolution$nrow,ncols = resolution$ncol),
        n = 1
      )
    })
  }
  brick(fieldlist)
}

tm <- Sys.time()
cat("Doing brick, time = ",format(tm),"\n")
# Do it for only 100, because it takes a lot of time
set.seed(123)
simstodo <- sample(ncol(posterior_samples$samps),100,replace = FALSE)
fieldsbrick <- simulate_spatial_fields(
  U = posterior_samples$samps[(ncol(BX) + 1):Wd,simstodo],
  theta = posterior_samples$theta[simstodo, 1:2],
  pointsdata = pointsdata,
  resolution = list(nrow = 400,ncol = 200)
)
cat("Finished brick, took ",format(difftime(Sys.time(),tm,units = 'secs')),"\n")
save(fieldsbrick,file = "spatial.Rdata")



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
simfieldsexceedence <- mean(fieldsbrick > log(1.2))

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




