repeating_smoothing_M_times <- function(M = 100, truth = "complicated", N = 1000 , baseline){
  result <- tibble()
  for(mm in 1:M){
    data <- Simulate_data_extreme(baseline = baseline, truth = truth, N = N)
    data <- abcoxph:::arrange_data(data)
    dat <- tibble(x = data$exposure, t = data$times, cens = data$censoring)
    dat$ranks <- rank(dat$t, ties.method = "min")
    ## setup smoothing part:
    a <- min(dat$x)
    b <- max(dat$x) # boundary
    n <- nrow(dat)
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
    P <- as(construct_penalty(dat$x,splineknots,p,m),'dgTMatrix')
    B <- as(construct_design(dat$x,splineknots,p,m),'dgTMatrix')
    D <- abcoxph:::create_diff_matrix(n)
    
    ### Setup TMB:
    tmbdat <- list(
      # Design matrix
      BX = B,
      # Penalty matrix
      P = P,
      # Differencing matrix
      D = D,
      # Log determinant of penalty matrix (without the sigma part)
      logPdet = as.numeric(determinant(P,logarithm = TRUE)$modulus),
      # Response
      ranks = as.integer(dat$ranks),
      cens = as.integer(dat$cens),
      # Prior params
      u = 2,
      alpha = 0.5
    )
    
    tmbparams <- list(
      W = rep(0,ncol(B)), # W = c(U); U = B-Spline coefficients
      theta = 0 # -2log(sigma)
    )
    ff <- TMB::MakeADFun(
      data = tmbdat,
      parameters = tmbparams,
      random = "W",
      DLL = "03_coxph_bspline",
      silent = TRUE
    )
    # Hessian not implemented for RE models
    ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
    # AGHQ
    quad <- aghq::marginal_laplace_tmb(ff,12,0)
    samps <- sample_marginal(quad,2000)
    W <- apply(samps$samps,1,mean)
    U <- W[1:ncol(P)]
    truefunc <- function(x) 1.5*(sin(0.8*x))
    plotx <- seq(a,b,by=0.01)
    sampled_lambda <- construct_design(plotx,splineknots,p,m) %*% samps$samps
    for (i in 1:ncol(sampled_lambda)) {
      sampled_lambda[,i] <- sampled_lambda[,i] - sampled_lambda[1,i]
    }
    est_lambda <- apply(sampled_lambda,1,mean)
    est_lambda <- est_lambda
    lambda_lower <- apply(sampled_lambda,1,quantile,probs = .025)
    lambda_upper <- apply(sampled_lambda,1,quantile,probs = .975)
    true_lambda <- truefunc(plotx) - truefunc(plotx)[1]
    rmse_aghq <- sqrt( mean( (est_lambda - true_lambda)^2 ) )
    mse_aghq <- mean( (est_lambda - true_lambda)^2 )
    
    covr_aghq <- mean(true_lambda <= lambda_upper & true_lambda >= lambda_lower)
    ### INLA:
    cnsA1 <- matrix(rep(0,50),nrow = 1)
    cnsA1[1] <- 1
    conse <- matrix(0, nrow = 1, ncol = 1)
    prior.prec <- list(prec = list(prior = "pc.prec",
                                   param = c(tmbdat$u, tmbdat$a)))
    formula <- inla.surv(times,censoring) ~ f(exposure_binned,model = 'rw2',constr = F, extraconstr = list(A=cnsA1,e=conse), hyper = prior.prec)
    Inlaresult <- inla(formula = formula, data = data, family = "coxph")
    fhat <- Inlaresult$summary.random$exposure_binned$mean
    fup <- Inlaresult$summary.random$exposure_binned$`0.975quant`
    flo <- Inlaresult$summary.random$exposure_binned$`0.025quant`
    fhat[1] = 0
    fup[1] = 0
    flo[1] = 0
    plotINLA <- data.frame(x = Inlaresult$summary.random$exposure_binned$ID, f = fhat, up = fup, lo = flo)
    plotINLA$true <- truefunc(plotINLA$x) - truefunc(plotINLA$x)[1]
    rmse_inla <- sqrt( mean( ((plotINLA$f) - (plotINLA$true))^2 ) )
    mse_inla <- mean( ((plotINLA$f) - (plotINLA$true))^2 )
    
    covr_inla <- mean((plotINLA$true) <= (plotINLA$up) & (plotINLA$true) >= (plotINLA$lo))
    ## mgcv:
    mgcvmod_bs <- mgcv::gam(
      t ~ s(x,bs='bs',m=c(p-1,m),k=length(splineknots)-p, pc = a),
      data = dat,
      family = cox.ph()
    )
    mgcv_bs_pred <- predict(mgcvmod_bs,newdata = data.frame(x = plotx),se.fit = TRUE)
    rmse_mgcv <- sqrt( mean( (((mgcv_bs_pred$fit) - ((truefunc(plotx)-truefunc(plotx)[1])))^2)[-1]) )
    mse_mgcv <- mean( (((mgcv_bs_pred$fit) - ((truefunc(plotx)-truefunc(plotx)[1])))^2)[-1])
    covr_mgcv <- mean((((truefunc(plotx)-truefunc(plotx)[1])) <= (mgcv_bs_pred$fit + 2*mgcv_bs_pred$se.fit) & ((truefunc(plotx)-truefunc(plotx)[1])) >= (mgcv_bs_pred$fit - 2*mgcv_bs_pred$se.fit))[-1])
    resultnew <- data.frame(
      method = c('AGHQ','MGCV','INLA'),
      rmse = c(rmse_aghq,rmse_mgcv,rmse_inla),
      coverage = c(covr_aghq,covr_mgcv,covr_inla),
      mse = c(mse_aghq, mse_mgcv, mse_inla)
    )
    result <- rbind(result,resultnew)
  }
  result
}