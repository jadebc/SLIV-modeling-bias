###################################################
# Shoo the Flu evaluation bias analysis 

# define priors for bias analysis that corrects
# flu hospitalization case counts for incomplete
# influenza testing

# output: prior distributions
###################################################

rm(list =ls())

# load libraries, source file paths
source(paste0(here::here(), "/bias-analysis/0-config.R"))

# define testing rate priors ---------------------------------------
# assume same by sex, race
# different correction priors for age, season
set.seed(123)
n = 10000

## age 0-4 ---------------------------------------
# non-kaiser
params_dis_0_4_nonk = find_beta_shape_params(mu = 0.635, sd = 0.025)
prior_dis_0_4_nonk = rbeta(n=n, 
                           shape1 = params_dis_0_4_nonk$a, 
                           shape2 = params_dis_0_4_nonk$b)
hist(prior_dis_0_4_nonk, breaks=100)
quantile(prior_dis_0_4_nonk, probs = c(0,0.5,1))

# kaiser  
params_dis_0_4_kaiser = find_beta_shape_params(mu = 0.78, sd = 0.07)
prior_dis_0_4_kaiser = rbeta(n=n, 
                             shape1 = params_dis_0_4_kaiser$a, 
                             shape2 = params_dis_0_4_kaiser$b)
hist(prior_dis_0_4_kaiser, breaks=100)
quantile(prior_dis_0_4_kaiser, probs = c(0,0.5,1))

# combined prior 
alpha_0_4 = params_dis_0_4_nonk$a * 0.7 + params_dis_0_4_kaiser$a * 0.3
beta_0_4 = params_dis_0_4_nonk$b * 0.7 + params_dis_0_4_kaiser$b * 0.3

prior_dis_0_4 = rbeta(n = n,
                      shape1 = alpha_0_4,
                      shape2 = beta_0_4)
hist(prior_dis_0_4, breaks=100)
quantile(prior_dis_0_4, probs = c(0,0.5,1))

## age 5-12 ---------------------------------------
# non-kaiser
params_dis_5_12_nonk = find_beta_shape_params(mu = 0.446, sd = 0.03)
prior_dis_5_12_nonk = rbeta(n=n, 
                            shape1 = params_dis_5_12_nonk$a, 
                            shape2 = params_dis_5_12_nonk$b)
hist(prior_dis_5_12_nonk, breaks=100)
quantile(prior_dis_5_12_nonk, probs = c(0,0.5,1))

# kaiser  
params_dis_5_12_kaiser = find_beta_shape_params(mu = 0.7, sd = 0.045)
prior_dis_5_12_kaiser = rbeta(n=n, 
                              shape1 = params_dis_5_12_kaiser$a, 
                              shape2 = params_dis_5_12_kaiser$b)
hist(prior_dis_5_12_kaiser, breaks=100)
quantile(prior_dis_5_12_kaiser, probs = c(0,0.5,1))

# combined prior 
alpha_5_12 = params_dis_5_12_nonk$a * 0.7 + params_dis_5_12_kaiser$a * 0.3
beta_5_12 = params_dis_5_12_nonk$b * 0.7 + params_dis_5_12_kaiser$b * 0.3

prior_dis_5_12 = rbeta(n = n,
                       shape1 = alpha_5_12,
                       shape2 = beta_5_12)
hist(prior_dis_5_12, breaks=100)
quantile(prior_dis_5_12, probs = c(0,0.5,1))


## age 13-17  ---------------------------------------
# non-kaiser
params_dis_13_17_nonk = find_beta_shape_params(mu = 0.537, sd = 0.025)
prior_dis_13_17_nonk = rbeta(n=n, 
                             shape1 = params_dis_13_17_nonk$a, 
                             shape2 = params_dis_13_17_nonk$b)
hist(prior_dis_13_17_nonk, breaks=100)
quantile(prior_dis_13_17_nonk, probs = c(0,0.5,1))

# kaiser  
params_dis_13_17_kaiser = find_beta_shape_params(mu = 0.712, sd = 0.045)
prior_dis_13_17_kaiser = rbeta(n=n, 
                               shape1 = params_dis_13_17_kaiser$a, 
                               shape2 = params_dis_13_17_kaiser$b)
hist(prior_dis_13_17_kaiser, breaks=100)
quantile(prior_dis_13_17_kaiser, probs = c(0,0.5,1))

# combined prior 
alpha_13_17 = params_dis_13_17_nonk$a * 0.7 + params_dis_13_17_kaiser$a * 0.3
beta_13_17 = params_dis_13_17_nonk$b * 0.7 + params_dis_13_17_kaiser$b * 0.3

prior_dis_13_17 = rbeta(n = n,
                        shape1 = alpha_13_17,
                        shape2 = beta_13_17)
hist(prior_dis_13_17, breaks=100)
quantile(prior_dis_13_17, probs = c(0,0.5,1))

## age 18-64  ---------------------------------------
# non-kaiser
params_dis_18_64_nonk = find_beta_shape_params(mu = 0.514, sd = 0.06)
prior_dis_18_64_nonk = rbeta(n=n, 
                             shape1 = params_dis_18_64_nonk$a, 
                             shape2 = params_dis_18_64_nonk$b)
hist(prior_dis_18_64_nonk, breaks=100)
quantile(prior_dis_18_64_nonk, probs = c(0,0.5,1))

# kaiser  
params_dis_18_64_kaiser = find_beta_shape_params(mu = 0.65, sd = 0.02)
prior_dis_18_64_kaiser = rbeta(n=n, 
                               shape1 = params_dis_18_64_kaiser$a, 
                               shape2 = params_dis_18_64_kaiser$b)
hist(prior_dis_18_64_kaiser, breaks=100)
quantile(prior_dis_18_64_kaiser, probs = c(0,0.5,1))

# combined prior 
alpha_18_64 = params_dis_18_64_nonk$a * 0.7 + params_dis_18_64_kaiser$a * 0.3
beta_18_64 = params_dis_18_64_nonk$b * 0.7 + params_dis_18_64_kaiser$b * 0.3

prior_dis_18_64 = rbeta(n = n,
                        shape1 = alpha_18_64,
                        shape2 = beta_18_64)
hist(prior_dis_18_64, breaks=100)
quantile(prior_dis_18_64, probs = c(0,0.5,1))


## age 65+  ---------------------------------------
# non-kaiser
params_dis_65_nonk = find_beta_shape_params(mu = 0.481, sd = 0.07)
prior_dis_65_nonk = rbeta(n=n, 
                          shape1 = params_dis_65_nonk$a, 
                          shape2 = params_dis_65_nonk$b)
hist(prior_dis_65_nonk, breaks=100)
quantile(prior_dis_65_nonk, probs = c(0,0.5,1))

# kaiser  
params_dis_65_kaiser = find_beta_shape_params(mu = 0.659, sd = 0.04)
prior_dis_65_kaiser = rbeta(n=n, 
                            shape1 = params_dis_65_kaiser$a, 
                            shape2 = params_dis_65_kaiser$b)
hist(prior_dis_65_kaiser, breaks=100)
quantile(prior_dis_65_kaiser, probs = c(0,0.5,1))

# combined prior 
alpha_65 = params_dis_65_nonk$a * 0.7 + params_dis_65_kaiser$a * 0.3
beta_65 = params_dis_65_nonk$b * 0.7 + params_dis_65_kaiser$b * 0.3

prior_dis_65 = rbeta(n = n,
                     shape1 = alpha_65,
                     shape2 = beta_65)
hist(prior_dis_65, breaks=100)
quantile(prior_dis_65, probs = c(0,0.5,1))

# save priors
prior_list <- list(
  prior_dis_0_4 = prior_dis_0_4,
  prior_dis_5_12 = prior_dis_5_12,
  prior_dis_13_17 = prior_dis_13_17,
  prior_dis_18_64 = prior_dis_18_64,
  prior_dis_65 = prior_dis_65
)
write.csv(prior_list, file = paste0(here::here(), "/results/priors_hosp.csv"))


# summarize priors
all_priors = bind_rows(
  quantile(prior_dis_0_4, probs = c(0,0.5,1)),
  quantile(prior_dis_5_12, probs = c(0,0.5,1)),
  quantile(prior_dis_13_17, probs = c(0,0.5,1)),
  quantile(prior_dis_18_64, probs = c(0,0.5,1)),
  quantile(prior_dis_65, probs = c(0,0.5,1))
) %>% mutate(age = c("0-4 years", "5-12 years", "13-17 years",
                     "18-64 years", "65+ years"),
             alpha = round(c(alpha_0_4, alpha_5_12, alpha_13_17,
                             alpha_18_64, alpha_65)),
             beta = round(c(beta_0_4, beta_5_12, beta_13_17,
                            beta_18_64, beta_65))) %>% 
  mutate(prior = paste0("Beta(alpha=", alpha, ", beta=", beta, ")")) %>% 
  dplyr::select(age, prior, everything())

write.csv(all_priors, file = paste0(here::here(), "/results/priors_hosp_table.csv"))
