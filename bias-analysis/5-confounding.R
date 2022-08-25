###################################################
# Shoo the Flu evaluation bias analysis 

# correct DID and relative reduction
# for time-dependent confounding

# probabilistic bias analysis

# output: relative reduction in DID corrected for confounding 
###################################################

rm(list =ls())

# load libraries, source file paths
source(paste0(here::here(), "/bias-analysis/0-config.R"))

# load data   ---------------------------------------
d_full <- readRDS(data_path)
d <- d_full$fluseasCDPH_2_5

# only keep years for 2017-18 analysis
d <- lapply(d, function(x) 
  x %>% filter(seas %in% c(1112, 1213, 1314, 1718)))

# load hospitalization priors
hosp_priors <- read.csv(paste0(here::here(), "/results/priors_hosp.csv"))


# correct hospitalization counts  -------------------------------------
corr_cases <- prob_correct_cases(data = d$eld,
                                 agecat = "eld",
                                 prior_dis_0_4 = hosp_priors$prior_dis_0_4,
                                 prior_dis_5_12 = hosp_priors$prior_dis_5_12,
                                 prior_dis_13_17 = hosp_priors$prior_dis_13_17,
                                 prior_dis_18_64 = hosp_priors$prior_dis_18_64,
                                 prior_dis_65 = hosp_priors$prior_dis_65) 


# define confounding priors -------------------------------------
# realistic scenario 1 -------------------------------------
set.seed(123)
RD_cd_pre_prior <- runif(n=10000, min=-0.0007 , max=0.0007)
RD_cd_post_prior <-runif(n=10000, min=-0.0053, max=0.0053)
PCE1_pre_prior <- rnorm(n=10000, mean=0.8, sd=0.01)
PCE0_pre_prior <- rnorm(n=10000, mean=0.8, sd=0.01)
PCE0_post_prior <- rnorm(n=10000, mean=0.9, sd=0.01)
PCE1_post_prior <- rnorm(n=10000, mean=0.9, sd=0.01)

# get quantiles of prior distributions
quantile(RD_cd_pre_prior, probs = c(0, 0.5, 1)) %>% formatC(format = "f") 
quantile(RD_cd_pre_prior, probs = c(0, 0.5, 1)) %>% formatC(format = "f") 
quantile(PCE1_pre_prior, probs = c(0, 0.5, 1)) %>% formatC(format = "f") 
quantile(PCE0_pre_prior, probs = c(0, 0.5, 1)) %>% formatC(format = "f") 
quantile(PCE0_post_prior, probs = c(0, 0.5, 1)) %>% formatC(format = "f") 
quantile(PCE1_post_prior, probs = c(0, 0.5, 1)) %>% formatC(format = "f") 

## run bias analysis -------------------------------------
set.seed(123)
results_eld <- replicate(
  n = 10000,
  prob_correct_hosp_confounding(
    d = d$eld,
    agecat = "eld",
    hosp_priors,
    RD_cd_pre_prior = RD_cd_pre_prior,
    PCE1_pre_prior = PCE1_pre_prior,
    PCE0_pre_prior = PCE0_pre_prior,
    RD_cd_post_prior = RD_cd_post_prior,
    PCE1_post_prior = PCE1_post_prior,
    PCE0_post_prior = PCE0_post_prior, correct_hosp=F
  ), simplify = F
) 

res = bind_rows(results_eld) %>% 
  mutate(z = ifelse(percent_red_DID_corrected <= 0 & percent_red_DID_corrected >= -0.02, 1, 0)) %>% 
  mutate(z = as.factor(z)) %>% 
  mutate(check_ousd_pre = 
           check_valid(corr_cases, time = "pre",
                       RD_cd = prior_RD_cd_pre, 
                       PCE1 = prior_PCE1_pre, 
                       PCE0 = prior_PCE0_pre)$ousd,
         check_wcc_pre = 
           check_valid(corr_cases, time = "pre",
                       RD_cd = prior_RD_cd_pre, 
                       PCE1 = prior_PCE1_pre, 
                       PCE0 = prior_PCE0_pre)$wcc,
         check_ousd_post = 
           check_valid(corr_cases, time = "post",
                       RD_cd = prior_RD_cd_post, 
                       PCE1 = prior_PCE1_post, 
                       PCE0 = prior_PCE0_post)$ousd,
         check_wcc_post = 
           check_valid(corr_cases, time = "post",
                       RD_cd = prior_RD_cd_post, 
                       PCE1 = prior_PCE1_post, 
                       PCE0 = prior_PCE0_post)$wcc) 

table(res$check_ousd_pre)
table(res$check_wcc_pre)
table(res$check_ousd_post)
table(res$check_wcc_post)

table(res$z) 

res = res %>% filter(check_ousd_pre=="pass" & 
                       check_wcc_pre=="pass" & 
                       check_ousd_post=="pass" & 
                       check_wcc_post=="pass")

table(res$z) 

# Bias-corrected estimate of intervention impact 
mean(res$DID_corrected)*100000
min(res$DID_corrected)*100000
max(res$DID_corrected)*100000

mean(res$percent_red_DID_corrected) 
min(res$percent_red_DID_corrected) 
max(res$percent_red_DID_corrected) 

# N simulations
nrow(res)

saveRDS(results_eld, paste0(here::here(), "/results/confounding-results-eld-scenario1.RDS"))
saveRDS(res, paste0(here::here(), "/results/confounding-results-eld-scenario1-processed.RDS"))

# extreme scenario 1 -------------------------------------
set.seed(123)
RD_cd_pre_prior <- runif(n=10000, min=-0.0007 , max=0.0007)
RD_cd_post_prior <-runif(n=10000, min=-0.0053, max=0.0053)
PCE1_pre_prior <- rbeta(n=10000,2,100)
PCE0_pre_prior <- rnorm(n=10000, mean=0.9, sd=0.01)
PCE0_post_prior <- rbeta(n=10000,2,100)
PCE1_post_prior <- rnorm(n=10000, mean=0.9, sd=0.01)

# get quantiles of prior distributions
quantile(RD_cd_pre_prior, probs = c(0, 0.5, 1)) %>% formatC(format = "f") 
quantile(RD_cd_pre_prior, probs = c(0, 0.5, 1)) %>% formatC(format = "f") 
quantile(PCE1_pre_prior, probs = c(0, 0.5, 1)) %>% formatC(format = "f") 
quantile(PCE0_pre_prior, probs = c(0, 0.5, 1)) %>% formatC(format = "f") 
quantile(PCE0_post_prior, probs = c(0, 0.5, 1)) %>% formatC(format = "f") 
quantile(PCE1_post_prior, probs = c(0, 0.5, 1)) %>% formatC(format = "f") 


## run bias analysis -------------------------------------
set.seed(123)
results_eld_ex1 <- replicate(
  n = 10000,
  prob_correct_hosp_confounding(
    d = d$eld,
    agecat = "eld",
    hosp_priors,
    RD_cd_pre_prior = RD_cd_pre_prior,
    PCE1_pre_prior = PCE1_pre_prior,
    PCE0_pre_prior = PCE0_pre_prior,
    RD_cd_post_prior = RD_cd_post_prior,
    PCE1_post_prior = PCE1_post_prior,
    PCE0_post_prior = PCE0_post_prior, 
    correct_hosp=F
  ), simplify = F
) 

res_ex1 = bind_rows(results_eld_ex1) %>% 
  mutate(z = ifelse(percent_red_DID_corrected <= 0 & percent_red_DID_corrected >= -0.02, 1, 0)) %>% 
  mutate(z = as.factor(z)) %>% 
  mutate(check_ousd_pre = 
           check_valid(corr_cases, time = "pre",
                       RD_cd = prior_RD_cd_pre, 
                       PCE1 = prior_PCE1_pre, 
                       PCE0 = prior_PCE0_pre)$ousd,
         check_wcc_pre = 
           check_valid(corr_cases, time = "pre",
                       RD_cd = prior_RD_cd_pre, 
                       PCE1 = prior_PCE1_pre, 
                       PCE0 = prior_PCE0_pre)$wcc,
         check_ousd_post = 
           check_valid(corr_cases, time = "post",
                       RD_cd = prior_RD_cd_post, 
                       PCE1 = prior_PCE1_post, 
                       PCE0 = prior_PCE0_post)$ousd,
         check_wcc_post = 
           check_valid(corr_cases, time = "post",
                       RD_cd = prior_RD_cd_post, 
                       PCE1 = prior_PCE1_post, 
                       PCE0 = prior_PCE0_post)$wcc) 

table(res_ex1$check_ousd_pre)
table(res_ex1$check_wcc_pre)
table(res_ex1$check_ousd_post)
table(res_ex1$check_wcc_post)

table(res_ex1$z) 

res_ex1 = res_ex1 %>% filter(check_ousd_pre=="pass" & 
                       check_wcc_pre=="pass" & 
                       check_ousd_post=="pass" & 
                       check_wcc_post=="pass")

table(res_ex1$z) 
prop.table(table(res_ex1$z) )

# Bias-corrected estimate of intervention impact 
mean(res_ex1$DID_corrected)*100000
min(res_ex1$DID_corrected)*100000
max(res_ex1$DID_corrected)*100000

# bias-corrected percent reduction 
mean(res_ex1$percent_red_DID_corrected) 
min(res_ex1$percent_red_DID_corrected) 
max(res_ex1$percent_red_DID_corrected) 

# N simulations
nrow(res_ex1)

## make plot -------------------------------------

# subset that matches model
plot_matched <- ggpairs(res_match, columns=c(14:15, 17:18, 13, 16), 
        lower = list(continuous = wrap("points", alpha = 0.2))) 

# subset that doesn't match model 
plot_unmatched <- ggpairs(res_unmatch, columns=c(14:15, 17:18, 13, 16),
        lower = list(continuous = wrap("points", alpha = 0.05))) 

ggsave(plot_matched, filename=paste0(here::here(), "/figures/ggpair_ex1_match.png"),
       width=16, height=8)
ggsave(plot_unmatched, filename=paste0(here::here(), "/figures/ggpair_ex1_unmatch.png"),
       width=16, height=8)

saveRDS(results_eld_ex1, paste0(here::here(), "/results/confounding-results-eld-scenario2.RDS"))
saveRDS(res_ex1, paste0(here::here(), "/results/confounding-results-eld-scenario2-processed.RDS"))

# extreme scenario 2 -------------------------------------
set.seed(123)
RD_cd_pre_prior <- runif(n=10000, min=-0.01 , max=0.01)
RD_cd_post_prior <-runif(n=10000, min=-0.01, max=0.01)
PCE1_pre_prior <- rnorm(n=10000, mean=0.8, sd=0.01)
PCE0_pre_prior <- rnorm(n=10000, mean=0.8, sd=0.01)
PCE0_post_prior <- rnorm(n=10000, mean=0.9, sd=0.01)
PCE1_post_prior <- rnorm(n=10000, mean=0.9, sd=0.01)

## run bias analysis -------------------------------------
set.seed(123)
results_eld_ex2 <- replicate(
  n = 10000,
  prob_correct_hosp_confounding(
    d = d$eld,
    agecat = "eld",
    hosp_priors,
    RD_cd_pre_prior = RD_cd_pre_prior,
    PCE1_pre_prior = PCE1_pre_prior,
    PCE0_pre_prior = PCE0_pre_prior,
    RD_cd_post_prior = RD_cd_post_prior,
    PCE1_post_prior = PCE1_post_prior,
    PCE0_post_prior = PCE0_post_prior, 
    correct_hosp=F
  ), simplify = F
) 

res_ex2 = bind_rows(results_eld_ex2) %>% 
  mutate(z = ifelse(percent_red_DID_corrected <= 0 & percent_red_DID_corrected >= -0.02, 1, 0)) %>% 
  mutate(z = as.factor(z)) %>% 
  mutate(check_ousd_pre = 
           check_valid(corr_cases, time = "pre",
                       RD_cd = prior_RD_cd_pre, 
                       PCE1 = prior_PCE1_pre, 
                       PCE0 = prior_PCE0_pre)$ousd,
         check_wcc_pre = 
           check_valid(corr_cases, time = "pre",
                       RD_cd = prior_RD_cd_pre, 
                       PCE1 = prior_PCE1_pre, 
                       PCE0 = prior_PCE0_pre)$wcc,
         check_ousd_post = 
           check_valid(corr_cases, time = "post",
                       RD_cd = prior_RD_cd_post, 
                       PCE1 = prior_PCE1_post, 
                       PCE0 = prior_PCE0_post)$ousd,
         check_wcc_post = 
           check_valid(corr_cases, time = "post",
                       RD_cd = prior_RD_cd_post, 
                       PCE1 = prior_PCE1_post, 
                       PCE0 = prior_PCE0_post)$wcc) 

table(res_ex2$check_ousd_pre)
table(res_ex2$check_wcc_pre)
table(res_ex2$check_ousd_post)
table(res_ex2$check_wcc_post)

table(res_ex2$z) 

res_ex2 = res_ex2 %>% filter(check_ousd_pre=="pass" & 
                               check_wcc_pre=="pass" & 
                               check_ousd_post=="pass" & 
                               check_wcc_post=="pass")

table(res_ex2$z) 
prop.table(table(res_ex2$z) )


# Bias-corrected estimate of intervention impact 
mean(res_ex2$DID_corrected)*100000
min(res_ex2$DID_corrected)*100000
max(res_ex2$DID_corrected)*100000

# bias-corrected percent reduction 
mean(res_ex2$percent_red_DID_corrected) 
min(res_ex2$percent_red_DID_corrected) 
max(res_ex2$percent_red_DID_corrected) 

# N simulations
nrow(res_ex2)

saveRDS(results_eld_ex2, paste0(here::here(), "/results/confounding-results-eld-scenario3.RDS"))
saveRDS(res_ex2, paste0(here::here(), "/results/confounding-results-eld-scenario3-processed.RDS"))



