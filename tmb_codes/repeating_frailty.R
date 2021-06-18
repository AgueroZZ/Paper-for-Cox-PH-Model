repeating_frailties_M_times <- function(M = 100, beta, N, K, sd){
  result <- tibble()
  for(mm in 1:M){
  data <- Simulate_grouped_data(N = N, bas = "constant", K = K, beta = beta, sdtheta = sd)
  data <- abcoxph:::arrange_data(data)
  dat <- tibble(x = data$x, t = data$times, cens = data$censoring, group = data$group)
  dat$ranks <- rank(dat$t, ties.method = "min")
  X <- as(as.matrix(dat$x),"dgTMatrix")
  B <- as(abcoxph:::create_blist_element(u = dat$group)$B,"dgTMatrix")
  D <- as(abcoxph:::create_diff_matrix(n), "dgTMatrix") ### n = K * N
  ### Setup TMB:
  tmbdat <- list(
    # Design matrix (random and fixed)
    B = as(B,"dgTMatrix"),
    X = as(X,"dgTMatrix"),
    # Differencing matrix
    D = as(D,"dgTMatrix"),
    # Response
    ranks = as.integer(dat$ranks),
    cens = as.integer(dat$cens),
    # Prior params
    u = 1,
    alpha = 0.5,
    betaprec = 0.001)
  tmbparams <- list(
    W = rep(0,ncol(B)+ncol(X)),
    theta = 0 # -2log(sigma)
  )
  ##### Fitting:
  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    random = "W",
    DLL = "00_coxph_frailty",
    silent = TRUE
  )
  # Hessian not implemented for RE models
  ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
  # AGHQ
  quad <- aghq::marginal_laplace_tmb(ff,3,0)
  samps <- sample_marginal(quad,1000)
  beta_est <- samps$samps[(K+1),]
  post_means <- apply(samps$samps, 1, mean)
  post_sds <- apply(samps$samps, 1, sd)
  post_up <- apply(samps$samps, 1, quantile, probs = 0.975)
  post_lo <- apply(samps$samps, 1, quantile, probs = 0.025)
  post_sum_aghq <- data.frame(rbind(post_means,post_sds,post_up,post_lo))
  rownames(post_sum_aghq) <- c("mean", "sd", "upper", "lower")
  colnames(post_sum_aghq) <- c(c(1:K),"beta")
  beta_cov_aghq <- ifelse(beta <= post_sum_aghq[,K+1][3] & beta >= post_sum_aghq[,K+1][4], 1, 0)
  beta_mse_aghq <- (post_sum_aghq[,K+1][1] - beta)^2
  ### INLA:
  prior.prec <- list(prec = list(prior = "pc.prec",
                                 param = c(tmbdat$u, tmbdat$a)))
  formula <- inla.surv(t,cens)~ x + f(group, model = "iid", hyper = prior.prec)
  Inlaresult <- inla(formula = formula, control.fixed = list(prec = 0.001), data = dat, family = "coxph")
  beta_cov_inla <- ifelse(beta <= Inlaresult$summary.fixed[2,]$'0.975quant' & beta >= Inlaresult$summary.fixed[2,]$'0.025quant', 1, 0)
  beta_mse_inla <- (Inlaresult$summary.fixed[2,]$mean - beta)^2
  frailty <- data %>% select(c(group,true)) %>% arrange(group) %>% unique(by = group)
  frailty$AGHQ <- as.numeric(post_sum_aghq[1,][-(K+1)])
  frailty$INLA <- as.numeric(Inlaresult$summary.random$group$mean)
  ### Random effects: MSE
  frailty_mse_aghq <- mean((frailty$true - frailty$AGHQ)^2)
  frailty_mse_inla <- mean((frailty$true - frailty$INLA)^2)
  ### Random effects: coverage
  frailty$AGHQ_up <- as.numeric(post_sum_aghq[3,][-(K+1)])
  frailty$INLA_up <- as.numeric(Inlaresult$summary.random$group$`0.975quant`)
  frailty$AGHQ_lo <- as.numeric(post_sum_aghq[4,][-(K+1)])
  frailty$INLA_lo <- as.numeric(Inlaresult$summary.random$group$`0.025quant`)
  ### AGHQ:
  frailty_cov_aghq <- mean((frailty$true >= frailty$AGHQ_lo & frailty$true <= frailty$AGHQ_up))
  ### INLA:
  frailty_cov_inla <- mean((frailty$true >= frailty$INLA_lo & frailty$true <= frailty$INLA_up))
  result <- rbind(result, c(beta_cov_aghq,beta_cov_inla,beta_mse_aghq,beta_mse_inla,frailty_cov_aghq,frailty_cov_inla,frailty_mse_aghq,frailty_mse_inla))
  }
  colnames(result) <- c("beta_cov_aghq","beta_cov_inla","beta_mse_aghq","beta_mse_inla","frailty_cov_aghq","frailty_cov_inla", "frailty_mse_aghq", "frailty_mse_inla")
  result
}