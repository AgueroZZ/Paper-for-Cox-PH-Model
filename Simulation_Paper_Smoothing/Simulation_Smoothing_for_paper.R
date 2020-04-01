source("~/Paper-for-Cox-PH-Model/codes/00-load-packages.R")
source("~/Paper-for-Cox-PH-Model/codes/01-functions.R")

cores <- 4
options(mc.cores = 4L)

set.seed(123)
tdom <- seq(0, 1000, by=0.001)
haz <- rep(0, length(tdom))
cut <- 40
for (i in 1:cut) {
  low <- as.numeric(quantile(tdom,(i-1)/cut))
  high <- as.numeric(quantile(tdom,(i)/cut))
  if(i %% 2 == 1){
    a <- runif(1,0,1)
    if(a > 0.3) haz[tdom<=high & tdom > low] <- 1/850
    else {
      c <- tdom[tdom<=high & tdom > low]
      haz[tdom<=high & tdom > low] <-(1/3000) *(c-min(c))
    }
  }
  if(i %% 2 == 0){
    a <- runif(1,0,1)
    if(a > 0.8){
      c <- tdom[tdom<=high & tdom > low]
      haz[tdom<=high & tdom > low] <- 1/100
    }
    else{
      haz[tdom<=high & tdom > low] <- sample(c(1/110,1/90),size = 1,prob = c(0.5,0.5))
    }
  }
}
plot(tdom, haz, type='l', xlab='Time', ylab='True baseline hazard')


true <- data_frame(time = tdom[tdom<=728.80400], hazard = haz[tdom<=728.80400])
ggplot(data = true) + aes(x = time, y = hazard) + geom_line() +theme_classic(base_size = 28) + labs(x = "", y = "")




N = 400
RW2BINS = 50
PARALLEL_EXECUTION = T

u <- runif(N)
x <- runif(N,min = -6, max = 6)
eta <- 1.5*sin(0.8*x) + 1.5 + rnorm(N,mean = 0 ,sd= exp(-.5*13))
truefunc <- function(x) 1.5*sin(0.8*x) + 1.5
tibble(x = c(-6,6)) %>%
  ggplot(aes(x = x)) +
  theme_light() +
  stat_function(fun = truefunc)


# Set up Model Specification:

failtimes <- c()
for (i in 1:N) {
  hazz <- haz * exp(eta[i])
  cumhaz <- cumsum(hazz*0.001)
  Surv <- exp(-cumhaz)
  failtimes[i] <- tdom[colSums(outer(Surv, u[i], `>`))]
}

# Set up Model Specification:

failtimes <- c()
for (i in 1:N) {
  hazz <- haz * exp(eta[i])
  cumhaz <- cumsum(hazz*0.001)
  Surv <- exp(-cumhaz)
  failtimes[i] <- tdom[colSums(outer(Surv, u[i], `>`))]
}



data <- data_frame(x=x,times = failtimes, entry = rep(0,length(length(u))),censoring = ifelse(failtimes>=1000,yes = 0, no=1))
for (i in 1:length(data$censoring)) {
  if(data$censoring[i]==1) data$censoring[i] <- rbinom(n=1,size=1,p=0.8)
}

data <- rename(data,exposure = x)
data <- data %>% as_tibble() %>%
  mutate(exposure_binned = bin_covariate(exposure,bins = RW2BINS,type = "quantile"))
data <- arrange_data(data)
data$ID <- 1:length(u)
Alist <- list()
Alist$exposure <- create_alist_element(data$exposure_binned)
model_data <- list(
  A = Alist,
  M = Alist %>% map("A") %>% map(ncol) %>% reduce(sum) - 1,
  n = length(unique(data$ID))
)
model_data$theta_logprior <- function(theta,prior_alpha = .5,prior_u = 2.5) {
  lambda <- -log(prior_alpha)/prior_u
  log(lambda/2) - lambda * exp(-theta/2) - theta/2
}

model_data$p <- 0
model_data$diffmat <- create_diff_matrix(model_data$n)
model_data$lambdainv <- create_full_dtcp_matrix(model_data$n)
model_data$Nd <- model_data$n - 1
model_data$Ne <- model_data$n
model_data$Wd <- model_data$M  + model_data$Nd
model_data$Wdf <- model_data$M  + model_data$Ne
model_data$times <- data$times
model_data$censoring <- data$censoring
model_data$entry <- data$entry
model_data$ID <- data$ID
model_data$diffmat <- create_diff_matrix(model_data$n)
model_data$lambdainv <- create_full_dtcp_matrix(model_data$n)
model_data$A$exposure$Ad <- model_data$diffmat %*% model_data$A$exposure$A


model_data$thetagrid <- mvQuad::createNIGrid(dim = 1,type = "GLe",level = 50)
mvQuad::rescale(model_data$thetagrid,domain = c(-3,3))
thetalist <- split(mvQuad::getNodes(model_data$thetagrid),rep(1:nrow(mvQuad::getNodes(model_data$thetagrid))))

# Random effect model specification data
model_data$modelspec <- model_data$A %>%
  purrr::map("model") %>%
  purrr::map2(.,names(.),~tibble(covariate = .y,model = .x)) %>%
  purrr::reduce(bind_rows)
model_data$vectorofcolumnstoremove <- round(RW2BINS/2)
cat("Finished creating model data!\n")


#Optimization:
control1 <- list(
  prec = 1e-08,
  stop.trust.radius = 1e-10,
  report.freq = 10,
  report.level = 4,
  start.trust.radius = 100,
  contract.threshold = .25,
  contract.factor = .5,
  expand.factor = 3,
  trust.iter = 2000,
  maxit = 3000,
  preconditioner = 0
)



tm <- proc.time()
sim1opt <- optimize_all_thetas_parallel(
  theta = model_data$thetagrid,
  model_data = model_data,
  optcontrol = control1,
  doparallel = PARALLEL_EXECUTION
)
rt <- proc.time() - tm



#Hyperparameter:Posterior
optresults_withlogpostpre <- add_log_posterior_values(sim1opt,model_data = model_data)
optresults_withlogpost <- normalize_optresults_logpost(optresults_withlogpostpre) 

margpost1 <- marginal_hyperparameter_posterior(1,optresults_withlogpost)
priorfunc <- function(x) exp(model_data$theta_logprior(x))
priorfuncsigma <- function(x) (2/x) * exp(model_data$theta_logprior(-2*log(x)))


thetapostplot1 <- margpost1$margpost %>%
  mutate(theta_post = exp(thetalogmargpost)) %>%
  ggplot(aes(x = theta)) +
  theme_classic() +
  geom_line(aes(y = theta_post),colour = "black",size = 0.5) +
  geom_line(aes(y = priorfunc(theta)),colour = "black",linetype = "dashed",size = 0.5) +
  # coord_cartesian(xlim = c(0,20)) +
  labs(y = "Density",x = "") +
  theme(text = element_text(size = PLOT_TEXT_SIZE))

sigmapostplot1 <- margpost1$margpost %>%
  mutate(sigma_post = exp(sigmalogmargpost)) %>%
  ggplot(aes(x = sigma)) +
  theme_classic() +
  geom_line(aes(y = sigma_post),colour = "black",size = 0.5) +
  geom_line(aes(y = priorfuncsigma(sigma)),colour = "black",linetype = "dashed",size = 0.5) +
  labs(x = "",y = "Density") +
  theme(text = element_text(size = PLOT_TEXT_SIZE))





#Final Comparison:

#Our Latent Field Computation:
tm2 <- proc.time()
margmeans_and_vars <- compute_marginal_means_and_variances(
  i = (model_data$Wd-model_data$M + 1):(model_data$Wd),
  model_results <- optresults_withlogpostpre,
  model_data <- model_data,
  lincomb = NULL,
  constrA = NULL
)
rt2 <- proc.time() - tm2
print(rt2)

margmeans_rw2 <- margmeans_and_vars$mean[1:(RW2BINS-1)]
margmeans_iid <- margmeans_and_vars$mean[RW2BINS:model_data$M]
margvars_rw2 <- margmeans_and_vars$variance[1:(RW2BINS-1)]
margsd <- sqrt(margvars_rw2)
vv <- model_data$vectorofcolumnstoremove

margmeanall <- c(
  margmeans_rw2[1:(vv-1)],
  0,
  margmeans_rw2[vv:length(margmeans_rw2)]
)
margsd <- c(
  margsd[1:(vv-1)],
  0,
  margsd[vv:length(margsd)]
)

#INLA:
cnsA1 <- matrix(rep(0,RW2BINS),nrow = 1)
cnsA1[model_data$vectorofcolumnstoremove] <- 1
conse <- matrix(0, nrow = 1, ncol = 1)

prior.prec <- list(prec = list(prior = "pc.prec",
                               param = c(2.5, 0.5)))

formula <- inla.surv(times,censoring) ~ f(exposure_binned,model = 'rw2',constr = F, extraconstr = list(A=cnsA1,e=conse), hyper = prior.prec)
Inlaresult <- inla(formula = formula, control.compute = list(dic=TRUE),control.inla = list(strategy = 'gaussian',int.strategy = 'grid', correct = FALSE),data = data, family = "coxph")
fhat <- Inlaresult$summary.random$exposure_binned$mean
fhat[model_data$vectorofcolumnstoremove] = 0
meanhere <- fhat
plotINLA <- data.frame(exposure = Inlaresult$summary.random$exposure_binned$ID)



#GAM:
data1 <- data
x <- sort(unique(model_data$A$exposure$u))
b <- gam(times~s(exposure,pc = plotINLA$exposure[vv]),
         family=cox.ph(),data=data1,weights=censoring)
y <- as.numeric(predict.gam(b,data_frame(exposure = x)))-as.numeric(predict.gam(b,data_frame(exposure=plotINLA$exposure[vv])))


#Plot:

#Plot:
simplot <- tibble(
  x = sort(unique(model_data$A$exposure$u)),
  mymean = margmeanall,
  mymeanlower = mymean - 2*margsd,
  mymeanupper = mymean + 2*margsd
) %>%
  ggplot(aes(x = x)) +
  theme_light() +
  geom_ribbon(aes(ymin = mymeanlower,ymax = mymeanupper),fill = "lightgrey",alpha = .5) +
  geom_line(aes(y = mymeanupper),colour = "black",linetype = "dotted") +
  geom_line(aes(y = mymeanlower),colour = "black",linetype = "dotted") +
  geom_line(aes(y = mymean),colour = 'black',linetype = 'dotdash') + 
  geom_line(aes(y = truefunc(x) - truefunc(x[vv])),colour = 'black',linetype = 'solid') + 
  xlab("") +
  ylab("") +
  theme_classic(base_size = 28)


new_compare <- simplot + geom_line(aes(y = meanhere),colour = "black",linetype = "dashed") 





datainla <- data_frame(x = Inlaresult$marginals.hyperpar$`Precision for exposure_binned`[,1], y = Inlaresult$marginals.hyperpar$`Precision for exposure_binned`[,2])
sigma <- 1/sqrt(datainla$x)
sigma_dens <- 2*(sigma^(-3))*datainla$y
datainla_sigma <- data_frame(sigma = sigma[-c(1:8)], sigma_dens = sigma_dens[-c(1:8)])

sigmapostplot1 <- margpost1$margpost %>%
  mutate(sigma_post = exp(sigmalogmargpost)) %>%
  ggplot(aes(x = sigma)) +
  theme_classic(base_size = 28) +
  geom_line(aes(y = sigma_post),colour = "black",linetype = "solid",size = 0.5) +
  labs(x = "",y = "") +
  geom_line(aes(y = priorfuncsigma(sigma)),colour = 'black',linetype = 'dashed',size = 0.5) + 
  geom_line(data = datainla_sigma, aes(y = sigma_dens, x = sigma),colour = 'black',linetype = 'dotdash',size = 0.5)


a <- brinla::bri.basehaz.plot(Inlaresult,plot = F)
a$basehaz <- exp(a$basehaz) * exp(Inlaresult$summary.fixed$mean)
a <- data_frame(time = a$time, baseline_hazard = a$basehaz)
ggplot(a, aes(time, baseline_hazard)) + geom_step() + theme_classic(base_size = 28) + labs(x = "", y = "")



plot(a$time,a$basehaz,type="n") 
segments(a$time[-length(a$time)],a$basehaz[-length(a$time)],a$time[-1],a$basehaz[-length(a$time)]) 





#Ploting:
ggsave(filename = "~/STA497/SmoothingSim_PosterTheta.pdf",plot = thetapostplot1)
ggsave(filename = "~/SmoothingSim_PosterSigma.pdf", plot = sigmapostplot1)

ggsave(filename = "~/INLA_baseline.pdf", plot = brinla::bri.basehaz.plot(Inlaresult))
 
ggsave(filename = "~/STA497/SmoothingSim_FinalPlot.pdf", plot = new_compare)
