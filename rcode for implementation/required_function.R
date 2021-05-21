#### Ploting function updated on Nov 28
library(ggplot2)
library(tidyverse)
library(latex2exp)

#### Plot of smoothing result with true function:
plot_smooth_withtrue <- function(proposed_model,truefun,TEXT_SIZE = 30){
  theme_set(theme_classic(base_size = 15))
  RW2BINS <- proposed_model$model_data$RW2BINS
  model_data <- proposed_model$model_data
  margmeanall <- proposed_model$marginal_latent$marginal_mean_smoothing
  margsd <- proposed_model$marginal_latent$marginal_sd_smoothing
  simplot <- dplyr::tibble(
    x = sort(unique(model_data$A[[1]]$u)),
    mymean = margmeanall,
    mymeanlower = mymean - 2*margsd,
    mymeanupper = mymean + 2*margsd
  ) %>%
    ggplot2::ggplot(ggplot2::aes(x = x)) +
    ggplot2::ylim(0,10) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = exp(mymeanlower),ymax = ifelse(exp(mymeanupper)>10,10,exp(mymeanupper))),fill = "orange",alpha = .2) +
    ggplot2::geom_line(ggplot2::aes(y = exp(mymean),color = "Fitted")) +
    geom_line(aes(y = exp(truefun(x) - truefun(x[round(RW2BINS/2)])),color = "True")) + labs(x = "u", y =  TeX('$exp(\\gamma(u)$)')) +
    scale_colour_manual(values=c("blue","black")) + ggplot2::theme(legend.position = c(.75, .92),legend.title = element_blank(),text = element_text(size = TEXT_SIZE),legend.background = element_blank())
  simplot
}


plot_smooth <- function(proposed_model,TEXT_SIZE = 30){
  theme_set(theme_classic(base_size = 15))
  RW2BINS <- proposed_model$model_data$RW2BINS
  model_data <- proposed_model$model_data
  margmeanall <- proposed_model$marginal_latent$marginal_mean_smoothing
  margsd <- proposed_model$marginal_latent$marginal_sd_smoothing
  simplot <- dplyr::tibble(
    x = sort(unique(model_data$A[[1]]$u)),
    mymean = margmeanall,
    mymeanlower = mymean - 2*margsd,
    mymeanupper = mymean + 2*margsd
  ) %>%
    ggplot2::ggplot(ggplot2::aes(x = x)) +
    ggplot2::ylim(0.5,1.5) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = exp(mymeanlower),ymax = ifelse(exp(mymeanupper)>10,10,exp(mymeanupper))),fill = "orange",alpha = .2) +
    ggplot2::geom_line(ggplot2::aes(y = exp(mymean),color = "Fitted")) +
    labs(x = "u", y =  TeX('$exp(\\gamma(u)$)')) + 
    scale_colour_manual(values=c("blue")) + ggplot2::theme(legend.position = c(.75, .92),legend.title = element_blank(),text = element_text(size = TEXT_SIZE),legend.background = element_blank())
  simplot
}



plot_hyper <- function(model, TEXT_SIZE = 30){
  theme_set(theme_classic(base_size = 15))
  getSig <- function(theta){
    sqrt(1/exp(theta))
  }
  getTheta <- function(sigma){
    log(1/(sigma^2))
  }
  translist <- list(totheta = getTheta, fromtheta = getSig)
  sig_pdf <- aghq::compute_pdf_and_cdf(model$thetaAGHQ$marginals[[1]],translist)
  priorfuncsigma <- function(x) (2/x) * exp(model$theta_logprior(-2*log(x)))
  psi0sig <- sig_pdf %>%
    ggplot(aes(x = transparam,y = pdf_transparam)) +
    theme_classic() +
    geom_line(aes(x = transparam,y = pdf_transparam,colour = 'Fitted')) +
    ggplot2::labs(x = TeX('$\\sigma'),y = "Density") +
    geom_line(aes(y = priorfuncsigma(transparam),colour = "Prior"), size = 0.5) +
    ggplot2::scale_colour_manual(values=c("blue","black")) +
    ggplot2::theme(legend.position = c(.75, .92),legend.title = element_blank(), text = element_text(size = TEXT_SIZE))
  psi0sig
}

plot_hyper_kid <- function(model,PL_sd, TEXT_SIZE = 30){
  theme_set(theme_classic(base_size = 15))
  getSig <- function(theta){
    sqrt(1/exp(theta))
  }
  getTheta <- function(sigma){
    log(1/(sigma^2))
  }
  translist <- list(totheta = getTheta, fromtheta = getSig)
  sig_pdf <- aghq::compute_pdf_and_cdf(model$thetaAGHQ$marginals[[1]],translist)
  priorfuncsigma <- function(x) (2/x) * exp(model$theta_logprior(-2*log(x)))
  psi0sig <- sig_pdf %>%
    ggplot(aes(x = transparam,y = pdf_transparam)) +
    theme_classic() +
    geom_line(aes(x = transparam,y = pdf_transparam,colour = 'Fitted')) +
    ggplot2::labs(x = TeX('$\\sigma'),y = "Density") +
    geom_line(aes(y = priorfuncsigma(transparam),colour = "Prior"), size = 0.5) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = PL_sd, colour = "PL Estimate")) + 
    ggplot2::scale_colour_manual(values=c("blue","red",'black')) +
    ggplot2::theme(legend.position = c(.75, .92),legend.title = element_blank(), text = element_text(size = TEXT_SIZE))
  psi0sig
}



### Simulation Example
Simulate_data <- function(N = 500, bas, truth, RW2BINS = 50){
  if(bas == "constant") {
    timelim <- 200
    tdom <- seq(0, timelim, by = 0.0001)
    haz <- rep(0, length(tdom))
    haz <- 0*tdom + 0.03
    true <- data.frame(time = tdom, hazard = haz)
    if(truth == "smooth") {
      u <- runif(N)
      x <- runif(N,min = 0, max = 6)
      truefunc <- function(x) log((x + 1)^2) - 1
      eta <- truefunc(x)
    }
    else{
      u <- runif(N)
      x <- runif(N,min = -5, max = 5)
      truefunc <- function(x) 1.5*(sin(0.8*x))
      eta <- truefunc(x)
    }
    failtimes <- c()
    for (i in 1:N) {
      hazz <- haz * exp(eta[i])
      cumhaz <- cumsum(hazz*0.0001)
      Surv <- exp(-cumhaz)
      failtimes[i] <- tdom[colSums(outer(Surv, u[i], `>`))]
    }
    data <- data_frame(x = x,times = failtimes, entry = rep(0,length(length(u))),censoring = ifelse(failtimes >= timelim,yes = 0, no=1))
    for (i in 1:length(data$censoring)) {
      if (data$censoring[i] == 1) {data$censoring[i] <- rbinom(n = 1,size = 1,prob = 0.9)}
    }
    data <- rename(data,exposure = x)
    data <- data %>% as_tibble() %>%
      mutate(exposure_binned = abcoxph:::bin_covariate(exposure,bins = RW2BINS,type = "quantile"))
  }
  else if (bas == "regular") {
    timelim <- 200
    tdom <- seq(0, timelim, by = 0.0001)
    haz <- rep(0, length(tdom))
    haz <- 0.03 * cos(0.15*tdom) + 0.05
    true <- data.frame(time = tdom, hazard = haz)
    if(truth == "smooth"){
      u <- runif(N)
      x <- runif(N,min = 0, max = 6)
      truefunc <- function(x) log((x + 1)^2) - 1
      eta <- truefunc(x)
    }
    else{
      u <- runif(N)
      x <- runif(N,min = -5, max = 5)
      truefunc <- function(x) 1.5*(sin(0.8*x))
      eta <- truefunc(x)
    }
    failtimes <- c()
    for (i in 1:N) {
      hazz <- haz * exp(eta[i])
      cumhaz <- cumsum(hazz*0.0001)
      Surv <- exp(-cumhaz)
      failtimes[i] <- tdom[colSums(outer(Surv, u[i], `>`))]
    }
    data <- data_frame(x = x,times = failtimes, entry = rep(0,length(length(u))),censoring = ifelse(failtimes >= timelim,yes = 0, no=1))
    for (i in 1:length(data$censoring)) {
      if (data$censoring[i] == 1) {data$censoring[i] <- rbinom(n = 1,size = 1,prob = 0.9)}
    }
    data <- rename(data,exposure = x)
    data <- data %>% as_tibble() %>%
      mutate(exposure_binned = abcoxph:::bin_covariate(exposure,bins = RW2BINS,type = "quantile"))
  }
  else{
    timelim <- 200
    tdom <- seq(0, timelim, by = 0.0001)
    haz <- rep(0, length(tdom))
    haz[1:round(length(haz)/3)] <- 0.01*sin(0.85*tdom[1:round(length(haz)/3)]) + 0.02
    haz[(round(length(haz)/3) + 1):(round(2*length(haz))/3)] <- 0.03 * cos(0.25*tdom[(round(length(haz)/3) + 1):(round(2*length(haz))/3)]) + 0.05
    haz[(round(2*length(haz)/3) + 1):length(haz)] <- 0.01*sin(0.55*tdom[(round(2*length(haz)/3) + 1):length(haz)]) + 0.025
    true <- data.frame(time = tdom, hazard = haz)
    if(truth == "smooth"){
      u <- runif(N)
      x <- runif(N,min = 0, max = 6)
      truefunc <- function(x) log((x + 1)^2) - 1
      eta <- truefunc(x)
    }
    else{
      u <- runif(N)
      x <- runif(N,min = -5, max = 5)
      truefunc <- function(x) 1.5*(sin(0.8*x))
      eta <- truefunc(x)
    }
    failtimes <- c()
    for (i in 1:N) {
      hazz <- haz * exp(eta[i])
      cumhaz <- cumsum(hazz*0.0001)
      Surv <- exp(-cumhaz)
      failtimes[i] <- tdom[colSums(outer(Surv, u[i], `>`))]
    }
    data <- data_frame(x = x,times = failtimes, entry = rep(0,length(length(u))),censoring = ifelse(failtimes >= timelim,yes = 0, no=1))
    for (i in 1:length(data$censoring)) {
      if (data$censoring[i] == 1) {data$censoring[i] <- rbinom(n = 1,size = 1,prob = 0.9)}
    }
    data <- rename(data,exposure = x)
    data <- data %>% as_tibble() %>%
      mutate(exposure_binned = abcoxph:::bin_covariate(exposure,bins = RW2BINS,type = "quantile"))
  }
  data
}


