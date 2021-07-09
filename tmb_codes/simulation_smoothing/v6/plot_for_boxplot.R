### load packages:
library(tidyverse)
TEXT_SIZE = 25

### load results:
load("/Users/ziangzhang/Documents/GitHub/Paper-for-Cox-PH-Model/tmb_codes/simulation_smoothing/v6/result_300_aggregations_complicated.rda")
load("/Users/ziangzhang/Documents/GitHub/Paper-for-Cox-PH-Model/tmb_codes/simulation_smoothing/v6/result_300_aggregations_stepwise.rda")
load("/Users/ziangzhang/Documents/GitHub/Paper-for-Cox-PH-Model/tmb_codes/simulation_smoothing/v6/result_300_aggregations_simple.rda")
# load("/home/ziang/Documents/Paper-for-Cox-PH-Model/tmb_codes/simulation_smoothing/v6/result_300_aggregations_complicated.rda")
# load("/home/ziang/Documents/Paper-for-Cox-PH-Model/tmb_codes/simulation_smoothing/v6/result_300_aggregations_stepwise.rda")
# load("/home/ziang/Documents/Paper-for-Cox-PH-Model/tmb_codes/simulation_smoothing/v6/result_300_aggregations_simple.rda")

### First plot: when true function is simple:
simple_result <- result3 %>% filter(method != "MGCV") %>% select(method,coverage,mse)
simple_result$method[simple_result$method == "AGHQ"] <- "Proposed"
simple_result %>% ggplot(aes(x = method, y = mse)) + geom_boxplot() + theme_classic(base_size = TEXT_SIZE)
simple_result <- simple_result %>% mutate(base = "simple")

### Second plot: when true function is stepwise:
stepwise_result <- result2 %>% filter(method != "MGCV") %>% select(method,coverage,mse)
stepwise_result$method[stepwise_result$method == "AGHQ"] <- "Proposed"
stepwise_result %>% ggplot(aes(x = method, y = mse)) + geom_boxplot() + theme_classic(base_size = TEXT_SIZE)
stepwise_result <- stepwise_result %>% mutate(base = "oscillating")


### Third plot: when true function is stepwise:
com_result <- result %>% filter(method != "MGCV") %>% select(method,coverage,mse)
com_result$method[com_result$method == "AGHQ"] <- "Proposed"
com_result %>% ggplot(aes(x = method, y = mse)) + geom_boxplot() + theme_classic(base_size = TEXT_SIZE)
com_result <- com_result %>% mutate(base = "complicated")

#### Side by side boxplot
result <- rbind(simple_result, stepwise_result, com_result)
result %>% ggplot(aes(x = base, y = mse, fill = method)) + geom_boxplot()  + 
  xlab("") + ylab("MSE") +
  theme_classic(base_size = TEXT_SIZE) +  theme(legend.position = "none")

#### barchart for coverage prob

result %>% group_by(base, method) %>% summarise(coverage = mean(coverage)) %>% 
  ggplot(aes(x = base, y = coverage, fill = method)) + geom_col(position = "dodge") +
  geom_hline(yintercept = 0.95, color = "red") +
  xlab("") + ylab("Coverage Prob") +
  theme_classic(base_size = TEXT_SIZE) +  theme(legend.position = "none")

