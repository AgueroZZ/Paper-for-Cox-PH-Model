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
# Set the resolution for the spatial interpolations.
reslist <- list(nrow = 50,ncol = 100)

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






