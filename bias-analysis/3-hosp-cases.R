###################################################
# Shoo the Flu evaluation bias analysis 

# correct hospitalization counts 
# for incomplete flu testing
# same priors in each site

# probabilistic bias analysis

# output: corrected case counts
###################################################

rm(list =ls())

# load libraries, source file paths
source(paste0(here::here(), "/bias-analysis/0-config.R"))

# load data   ---------------------------------------
d_full <- readRDS(data_path)
d <- d_full$fluseasCDPH_2_5

eld <- d$eld

# read in testing rate priors ---------------------------------------
all_priors <- read.csv(paste0(here::here(), "/results/priors_hosp.csv"))

# run probabilistic bias analysis ---------------------------------------
set.seed(123)
results_list = replicate(
  n = 10000, 
  prob_correct_cases(data = eld,
                     agecat= "eld",
                   prior_dis_0_4 = all_priors$prior_dis_0_4,
                   prior_dis_5_12 = all_priors$prior_dis_5_12,
                   prior_dis_13_17 = all_priors$prior_dis_13_17,
                   prior_dis_18_64 = all_priors$prior_dis_18_64,
                   prior_dis_65 = all_priors$prior_dis_65),
  simplify= F
)


results_df = bind_rows(results_list)

results_did = lapply(results_list, calculate_did)
results_did_vec = unlist(results_did)

# summarize results ---------------------------------------
mean(results_did_vec)
min(results_did_vec)
max(results_did_vec)

results_relred = lapply(results_list, calculate_relred)
results_relred_vec = unlist(results_relred)

mean(results_relred_vec)
min(results_relred_vec)
max(results_relred_vec)

hist(results_relred_vec)

saveRDS(results_df, paste0(here::here(), "/results/hosp-corrected-cases.RDS"))




