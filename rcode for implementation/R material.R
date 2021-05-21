### To Install AGHQ: devtools::install_github("awstringer1/aghq")
### Before you install Abcoxph package, make sure all the other dependancy packages have been installed.
### To Install Abcoxph: devtools::install_github("AgueroZZ/Abcoxph")


### Source the needed R package:
source("~/desktop/load_all.R")


### Source the required function for plotting and simulating of data:
source("~/desktop/required_function.R")


## constant baseline and complicated truth
set.seed(4321)
data <- Simulate_data(bas = "constant", truth = "complicated")
prior <- prior_setup_AGHQ(theta_number = 3, startingvals = 1,theta_prior_a = .5,theta_prior_u = 2.5)
system.time(model3 <- abcoxph_AGHQ_setup(times~s(exposure), cens = "censoring", data = data, prior_control = prior, RW2BINS = 50))
proposed_model3 <- abcox_fit(model3, PARALLEL_EXECUTION = F)
plot_smooth_withtrue(proposed_model3,truefun = function(x) {1.5*(sin(0.8*x))})
plot_hyper(model3)

## Complicated baseline and complicated truth
set.seed(321)
data <- Simulate_data(bas = "complicated", truth = "complicated")
prior <- prior_setup_AGHQ(theta_number = 3, startingvals = 1,theta_prior_a = .5,theta_prior_u = 2.5)
system.time(model1 <- abcoxph_AGHQ_setup(times~s(exposure), cens = "censoring", data = data, prior_control = prior, RW2BINS = 50))
proposed_model1 <- abcox_fit(model1, PARALLEL_EXECUTION = F)
plot_smooth_withtrue(proposed_model1,truefun = function(x) {1.5*(sin(0.8*x))})
plot_hyper(model1)

## regular baseline and complicated truth
set.seed(21)
data <- Simulate_data(bas = "regular", truth = "complicated")
prior <- prior_setup_AGHQ(theta_number = 3, startingvals = 1,theta_prior_a = .5,theta_prior_u = 2.5)
system.time(model2 <- abcoxph_AGHQ_setup(times~s(exposure), cens = "censoring", data = data, prior_control = prior, RW2BINS = 50))
proposed_model2 <- abcox_fit(model2, PARALLEL_EXECUTION = F)
plot_smooth_withtrue(proposed_model2,truefun = function(x) {1.5*(sin(0.8*x))})
plot_hyper(model2)





## Complicated baseline and smooth truth
set.seed(123)
data <- Simulate_data(bas = "complicated", truth = "smooth")
prior <- prior_setup_AGHQ(theta_number = 3, startingvals = 1,theta_prior_a = .5,theta_prior_u = 2.5)
system.time(model4 <- abcoxph_AGHQ_setup(times~s(exposure), cens = "censoring", data = data, prior_control = prior, RW2BINS = 50))
proposed_model4 <- abcox_fit(model4, PARALLEL_EXECUTION = F)
plot_smooth_withtrue(proposed_model4,truefun = function(x) {log((x + 1)^2) - 1})
plot_hyper(model4)

## regular baseline and smooth truth
set.seed(1234)
data <- Simulate_data(bas = "regular", truth = "smooth")
prior <- prior_setup_AGHQ(theta_number = 3, startingvals = 1,theta_prior_a = .5,theta_prior_u = 2.5)
system.time(model5 <- abcoxph_AGHQ_setup(times~s(exposure), cens = "censoring", data = data, prior_control = prior, RW2BINS = 50))
proposed_model5 <- abcox_fit(model5, PARALLEL_EXECUTION = F)
plot_smooth_withtrue(proposed_model5,truefun = function(x) {log((x + 1)^2) - 1})
plot_hyper(model5)

## constant baseline and smooth truth
set.seed(12345)
data <- Simulate_data(bas = "constant", truth = "smooth")
prior <- prior_setup_AGHQ(theta_number = 3, startingvals = 1,theta_prior_a = .5,theta_prior_u = 2.5)
system.time(model6 <- abcoxph_AGHQ_setup(times~s(exposure), cens = "censoring", data = data, prior_control = prior, RW2BINS = 50))
proposed_model6 <- abcox_fit(model6, PARALLEL_EXECUTION = F)
plot_smooth_withtrue(proposed_model6,truefun = function(x) {log((x + 1)^2) - 1})
plot_hyper(model6)










### Leukemia Example:
data <- INLA::Leuk
prior <- prior_setup_AGHQ(theta_number = 3, startingvals = 1,theta_prior_a = .5,theta_prior_u = 2)
system.time(model <- abcoxph_AGHQ_setup(time~age + sex + wbc + s(tpi), cens = "cens", data = data, prior_control = prior, RW2BINS = 50))
model_fit <- abcox_fit(model,PARALLEL_EXECUTION = F)
plot_smooth(model_fit)
plot_hyper(model)












### Kidney Example:
data <- survival::kidney
data <- data.frame(time = data$time, status = data$status, age = data$age,sex = data$sex,GN = ifelse(data$disease == "GN",1,0),AN = ifelse(data$disease == "AN",1,0),PKD = ifelse(data$disease == "PKD",1,0), id = data$id)
prior <- prior_setup_AGHQ(theta_number = 3, startingvals = 1,theta_prior_a = .5,theta_prior_u = 2)
system.time(model <- abcoxph_AGHQ_setup(time~age + sex + GN + AN + PKD + id(id), cens = "status",data = data, prior_control = prior, correction = T))
model_fit2 <- abcox_fit(model,PARALLEL_EXECUTION = F)
model_fit2$marginal_latent

coxphfit <- coxph(Surv(time, status) ~ age + sex + GN + AN + PKD + frailty(id,dist = "gauss", sparse = F), data = data,ties = "breslow")
coxphfit$coefficients
summary(coxphfit)

plot_hyper_kid(model,PL_sd = 0.695)



