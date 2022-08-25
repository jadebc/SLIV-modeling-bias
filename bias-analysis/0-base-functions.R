



##############################################
# Documentation:  get_did_rates
# Usage:          get_did_rates(data)
# Description:    Calculate a and b shape parameters for a beta
#                 distribution from mean and standard deviation
#
# Args/Options:
# data:           data frame with flu hospitalization counts and N by season

# Returns: a data frame with crude difference-in-differences
# Output:  none
##############################################
get_did_rates <- function(data){
  post = data %>% filter(seas==1718) %>% 
    group_by(dist) %>% 
    summarise(cases = sum(flucases),
              cases_corr = sum(flucases_corr),
              N = sum(N)) %>% 
    mutate(pre = "post")
  
  
  pre = data %>% filter(seas>=1112 & seas<=1314) %>% 
    group_by(dist, seas) %>% 
    summarise(cases = sum(flucases),
              cases_corr = sum(flucases_corr),
              N = sum(N)) %>% 
    group_by(dist) %>% 
    summarise(cases = sum(cases),
              cases_corr = sum(cases_corr),
              N= sum(N)) %>% 
    mutate(pre = "pre")
  
  all=bind_rows(pre, post) %>% 
    mutate(rate = cases/N*10000,
           rate_corr = cases_corr/N*10000)
  
  return(all)
}




##############################################
# Documentation:  find_beta_shape_params
# Usage:          find_beta_shape_params(mu, sd)
# Description:    Calculate a and b shape parameters for a beta
#                 distribution from mean and standard deviation
#
# Args/Options:
# mu:             mean of distribution
# sd:             standard deviation of distribution

# Returns: a list with elements a and b that include shape parameters
# Output:  none
##############################################
find_beta_shape_params = function(mu, sd){
  var = sd ^ 2
  a = ((1 - mu) / var - 1 / mu) * mu ^ 2
  b = a * (1 / mu - 1)
  return (list(a = a, b = b))
}




##############################################
# Documentation:  correct_cases
# Usage:          correct_cases(cases, test_rate)
# Description:    correct hospitalization rates  
#                 true cases = ((cases *(1-test rate))/ test rate) + cases
#                 function to correct hospitalization count for incomplete testing
#
# Args/Options:
# cases:          observed case count
# test_rate:      probability of testing for influenza

# Returns: case count corrected for incomplete testing
# Output:  none
##############################################
correct_cases <- function(cases, test_rate){
  return(round(((cases *(1-test_rate))/ test_rate) + cases,0))
}


##############################################
# Documentation:  prob_correct_cases
# Usage:          prob_correct_cases(data, agecat, prior_dis_0_4, 
#                                    prior_dis_5_12, prior_dis_13_17, 
#                                    prior_dis_18_64, prior_dis_65)
# Description:    probabilistically correct case count
#
# Args/Options:
# data:           data frame that contains uncorrected case counts
# agecat:         "eld" for 65+ years, "all" for all ages, "nonelem" for non-elementary aged
# prior_dis_0_4:  prior distribution for influenza testing rate in 0-4 year olds
# prior_dis_5_12:  prior distribution for influenza testing rate in 5-12 year olds
# prior_dis_13_17:  prior distribution for influenza testing rate in 13-17 year olds
# prior_dis_18_64:  prior distribution for influenza testing rate in 18-64 year olds
# prior_dis_65:  prior distribution for influenza testing rate in 65+ year olds

# Returns: case count corrected for incomplete testing
# Output:  none
##############################################
prob_correct_cases <- function(data, agecat, prior_dis_0_4, 
                               prior_dis_5_12, 
                               prior_dis_13_17, 
                               prior_dis_18_64, 
                               prior_dis_65){
  
  draw_prior_dis_0_4 = sample(prior_dis_0_4, size = 1)
  draw_prior_dis_5_12 = sample(prior_dis_5_12, size = 1)
  draw_prior_dis_13_17 = sample(prior_dis_13_17, size = 1)
  draw_prior_dis_18_64 = sample(prior_dis_18_64, size = 1)
  draw_prior_dis_65 = sample(prior_dis_65, size = 1)
  
  if(agecat == "all"){
    data <- data %>% mutate(flucases_corr = case_when(
      agecat %in% c("Under 5") ~ correct_cases(cases = flucases, 
                                               test_rate = draw_prior_dis_0_4),
      agecat %in% c("5 to 9", "10 to 14") ~ correct_cases(cases = flucases, 
                                                          test_rate = draw_prior_dis_5_12),
      agecat %in% c("15 to 17") ~ correct_cases(cases = flucases, 
                                                test_rate = draw_prior_dis_5_12),
      agecat %in% c("18 and 19", "20", "21", "22 to 24",
                    "25 to 29", "30 to 34", "35 to 39",
                    "40 to 44", "45 to 49", "50 to 54",
                    "55 to 59", "60 and 61", "62 to 64") ~ correct_cases(cases = flucases, 
                                                                         test_rate = draw_prior_dis_18_64),
      agecat %in% c("65 and 66", "67 to 69",
                    "70 to 74", "75 to 79",
                    "80 to 84", "85 and over") ~ correct_cases(cases = flucases, 
                                                               test_rate = draw_prior_dis_65))
    )
    
  }
  
  if(agecat == "nonelem"){
    data <- data %>% mutate(flucases_corr = case_when(
      agecat %in% c("Under 5") ~ correct_cases(cases = flucases, 
                                               test_rate = draw_prior_dis_0_4),
      agecat %in% c("15 to 17") ~ correct_cases(cases = flucases, 
                                                test_rate = draw_prior_dis_5_12),
      agecat %in% c("18 and 19", "20", "21", "22 to 24",
                    "25 to 29", "30 to 34", "35 to 39",
                    "40 to 44", "45 to 49", "50 to 54",
                    "55 to 59", "60 and 61", "62 to 64") ~ correct_cases(cases = flucases, 
                                                                         test_rate = draw_prior_dis_18_64),
      agecat %in% c("65 and 66", "67 to 69",
                    "70 to 74", "75 to 79",
                    "80 to 84", "85 and over") ~ correct_cases(cases = flucases, 
                                                               test_rate = draw_prior_dis_65))
    )
  }
  
  if(agecat == "eld"){
    data <- data %>% mutate(flucases_corr = 
                              correct_cases(cases = flucases, 
                                            test_rate = draw_prior_dis_65))
  }
  
  
  rates = get_did_rates(data)
  rates_df = rates %>% 
    mutate(analysis = agecat) %>% 
    dplyr::select(analysis, dist, pre, cases, cases_corr, N) 
  
  out = rates_df %>% 
    mutate(prior_dis_0_4 = draw_prior_dis_0_4,
           prior_dis_5_12 = draw_prior_dis_5_12,
           prior_dis_13_17 = draw_prior_dis_13_17,
           prior_dis_18_64 = draw_prior_dis_18_64,
           prior_dis_65 = draw_prior_dis_65)
  
  return(out)
}


##############################################
# Documentation:  prob_correct_hosp_confounding
# Usage:          prob_correct_hosp_confounding(d, 
#                 agecat, hosp_priors, RD_cd_pre_prior, 
#                 PCE1_pre_prior, PCE0_pre_prior,
#                 RD_cd_post_prior, PCE1_post_prior, 
#                 PCE0_post_prior, correct_hosp)
# Description:    probabilistically correct case count
#                 for under testing and for confounding 
#
# Args/Options:
# d:                data frame that contains uncorrected case counts
# agecat:           "eld" for 65+ years, "all" for all ages, "nonelem" for non-elementary aged
# hosp_priors       data frame that contains flu testing rate priors by age as columns 
# RD_cd_pre_prior:  prior distribution for risk difference for confounder & disease (pre-intervention)
# PCE1_pre_prior:   prior distribution for probability of confounder in intervention group (pre-intervention)
# PCE0_pre_prior:   prior distribution for probability of confounder in comparison group (pre-intervention)
# RD_cd_post_prior: prior distribution for risk difference for confounder & disease (post-intervention)
# PCE1_pre_prior:   prior distribution for probability of confounder in intervention group (post-intervention)
# PCE1_pre_prior:   prior distribution for probability of confounder in intervention group (post-intervention)
# correct_hosp:     boolean (T/F) to indicate whether to correct for flu testing rates among hospitalized  

# Returns: case count corrected for incomplete testing and time-dependent confounding
# Output:  none
##############################################
prob_correct_hosp_confounding <- function(d, 
                                          agecat,
                                          hosp_priors,
                                          RD_cd_pre_prior, PCE1_pre_prior, PCE0_pre_prior,
                                          RD_cd_post_prior, PCE1_post_prior, PCE0_post_prior,
                                          correct_hosp=T){
  
  # correct hospitalization counts
  corr_cases <- prob_correct_cases(data = d,
                                   agecat = agecat,
                                   prior_dis_0_4 = hosp_priors$prior_dis_0_4,
                                   prior_dis_5_12 = hosp_priors$prior_dis_5_12,
                                   prior_dis_13_17 = hosp_priors$prior_dis_13_17,
                                   prior_dis_18_64 = hosp_priors$prior_dis_18_64,
                                   prior_dis_65 = hosp_priors$prior_dis_65) 
  
  # correct for confounding 
  RD_cd_pre_draw = sample(RD_cd_pre_prior, size = 1)
  PCE1_pre_draw = sample(PCE1_pre_prior, size = 1)
  PCE0_pre_draw = sample(PCE0_pre_prior, size = 1)
  RD_cd_post_draw = sample(RD_cd_post_prior, size = 1)
  PCE1_post_draw = sample(PCE1_post_prior, size = 1)
  PCE0_post_draw = sample(PCE0_post_prior, size = 1)
  
  DID_corr = correct_confounding(
    corr_cases = corr_cases,
    RD_cd_pre = RD_cd_pre_draw,
    PCE1_pre = PCE1_pre_draw,
    PCE0_pre = PCE0_pre_draw,
    RD_cd_post = RD_cd_post_draw,
    PCE1_post = PCE1_post_draw,
    PCE0_post = PCE0_post_draw, 
    correct_hosp=correct_hosp
  )
  
  priors = data.frame(
    prior_hosp_dis_0_4 = corr_cases$prior_dis_0_4 %>% unique(),
    prior_hosp_dis_5_12 = corr_cases$prior_dis_5_12 %>% unique(),
    prior_hosp_dis_13_17 = corr_cases$prior_dis_13_17 %>% unique(),
    prior_hosp_dis_18_64 = corr_cases$prior_dis_18_64 %>% unique(),
    prior_hosp_dis_65 = corr_cases$prior_dis_65 %>% unique()
  )
  
  out = bind_cols(DID_corr, priors) %>% 
    rename(
      prior_RD_cd_pre = RD_cd_pre,
      prior_PCE1_pre = PCE1_pre,
      prior_PCE0_pre = PCE0_pre,
      prior_RD_cd_post = RD_cd_post,
      prior_PCE1_post = PCE1_post,
      prior_PCE0_post = PCE0_post)
  
  return(out)
  
}


##############################################
# Documentation:  check_valid
# Usage:          check_valid(corr_cases, time,
#                 RD_cd, PCE1, PCE0)
# Description:    check whether a set of corrected case count yields
#                 non-negative cell-values 
#
# Args/Options:
# corr_cases:     data frame that contains corrected case counts
# time:           "pre" for pre-intervention, "post" for post-intervention
# RD_cd:  prior distribution for risk difference for confounder & disease
# PCE1:   prior distribution for probability of confounder in intervention group
# PCE0:   prior distribution for probability of confounder in comparison group

# Returns: vector indicating whether any results were invalid 
# Output:  none
##############################################
check_valid <- function(corr_cases, time,
                        RD_cd, PCE1, PCE0){
  
  a  = corr_cases$cases_corr[corr_cases$dist=="OUSD" & corr_cases$pre==time]
  b  = corr_cases$cases_corr[corr_cases$dist=="WCCUSD" & corr_cases$pre==time]
  m  = corr_cases$N[corr_cases$dist=="OUSD" & corr_cases$pre==time]
  n  = corr_cases$N[corr_cases$dist=="WCCUSD" & corr_cases$pre==time]
  
  a_crude  = corr_cases$cases[corr_cases$dist=="OUSD" & corr_cases$pre==time]
  b_crude  = corr_cases$cases[corr_cases$dist=="WCCUSD" & corr_cases$pre==time]
  
  # ---- pre intervention ----
  M1 = PCE1 * m
  N1 = PCE0 * n
  M0 = m - M1
  N0 = n - N1
  
  A1 = (M1*(m - M1)*RD_cd + M1*a)/m
  B1 = (N1*(n - N1)*RD_cd + N1*b)/n
  A0 = a - A1
  B0 = b - B1
  
  any_zero_cases_ousd = ifelse(
    A1<0 | A0<0 ,
    "fail", "pass"
  )
  
  any_zero_cases_wcc = ifelse(
    B1<0 | B0<0,
    "fail", "pass"
  )
  
  return(list(ousd = any_zero_cases_ousd,
              wcc = any_zero_cases_wcc))
  
}


##############################################
# Documentation:  get_did_rates_dist
# Usage:          get_did_rates_dist(data, ousd_prior, wcc_prior)
# Description:    check whether a set of corrected case count yields
#                 non-negative cell-values 
#
# Args/Options:
# data:           data frame that contains uncorrected case counts
# ousd_prior:     probability of testing for influenza in intervention group
# wcc_prior:      probability of testing for influenza in comparison group

# Returns: corrected relative reduction accounting for pre-intervention differences
# Output:  none
##############################################
get_did_rates_dist = function(data, ousd_prior, wcc_prior){
  int <- d$eld %>% mutate(flucases_corr = 
                            correct_cases(cases = flucases, 
                                          test_rate = ousd_prior)) %>% 
    filter(dist=="OUSD")
  
  comp <- d$eld %>% mutate(flucases_corr = 
                             correct_cases(cases = flucases, 
                                           test_rate = wcc_prior)) %>% 
    filter(dist=="WCCUSD")
  
  post_int = int %>% filter(seas==1718) %>% 
    summarise(cases = sum(flucases),
              cases_corr = sum(flucases_corr),
              N = sum(N)) %>% 
    mutate(pre = "post") %>% 
    mutate(dist= "OUSD")
  
  post_comp = comp %>% filter(seas==1718) %>% 
    summarise(cases = sum(flucases),
              cases_corr = sum(flucases_corr),
              N = sum(N)) %>% 
    mutate(pre = "post") %>% 
    mutate(dist= "WCCUSD")
  
  
  pre_int = int %>% filter(seas>=1112 & seas<=1314) %>% 
    group_by(seas) %>% 
    summarise(cases = sum(flucases),
              cases_corr = sum(flucases_corr),
              N = sum(N)) %>% 
    summarise(cases = sum(cases),
              cases_corr = sum(cases_corr),
              N= sum(N)) %>% 
    mutate(pre = "pre",
           dist="OUSD")
  
  pre_comp = comp %>% filter(seas>=1112 & seas<=1314) %>% 
    group_by(seas) %>% 
    summarise(cases = sum(flucases),
              cases_corr = sum(flucases_corr),
              N = sum(N)) %>% 
    summarise(cases = sum(cases),
              cases_corr = sum(cases_corr),
              N= sum(N)) %>% 
    mutate(pre = "pre",
           dist="WCCUSD")
  
  all=bind_rows(pre_int, post_int, pre_comp, post_comp) %>% 
    mutate(rate = cases/N*10000,
           rate_corr = cases_corr/N*10000)
  
  # corrected did
  did = (all$rate_corr[all$dist=="OUSD" & all$pre=="post"] -
           all$rate_corr[all$dist=="OUSD" & all$pre=="pre"]) -
    (all$rate_corr[all$dist=="WCCUSD" & all$pre=="post"] -
       all$rate_corr[all$dist=="WCCUSD" & all$pre=="pre"])
  
  # corrected RD in control (post-pre)
  RD_cont = (all$rate_corr[all$dist=="WCCUSD" & all$pre=="post"] -
               all$rate_corr[all$dist=="WCCUSD" & all$pre=="pre"]) 
  
  # corrected relative reduction 
  # (RDMH, post - RDMH, pre)/ RDMH, comparison. 
  relred = did / RD_cont
  
  return(relred)
}



##############################################
# Documentation:  calculate_did
# Usage:          calculate_did(rates)
# Description:    calculate difference-in-difference from 
#                 corrected cases counts 
#
# Args/Options:
# data:           data frame that contains corrected case counts

# Returns: corrected difference-in-difference
# Output:  none
##############################################
calculate_did <- function(rates){
  
  rates <- rates %>% mutate(
    rate_corr = cases_corr/N,
    rate = cases/N
  )
  # corrected did
  ((rates$rate_corr[rates$dist=="OUSD" & rates$pre=="post"] -
      rates$rate_corr[rates$dist=="OUSD" & rates$pre=="pre"]) -
      (rates$rate_corr[rates$dist=="WCCUSD" & rates$pre=="post"] -
         rates$rate_corr[rates$dist=="WCCUSD" & rates$pre=="pre"]))*100000
}

##############################################
# Documentation:  calculate_relred
# Usage:          calculate_relred(rates)
# Description:    calculate relative reduction difference-in-difference from 
#                 corrected cases counts 
#
# Args/Options:
# data:           data frame that contains corrected case counts

# Returns: corrected relative reduction account for pre-intervention differences
# Output:  none
##############################################
calculate_relred <- function(rates){
  
  rates <- rates %>% mutate(
    rate_corr = cases_corr/N,
    rate = cases/N
  )
  # corrected did
  did = ((rates$rate_corr[rates$dist=="OUSD" & rates$pre=="post"] -
            rates$rate_corr[rates$dist=="OUSD" & rates$pre=="pre"]) -
           (rates$rate_corr[rates$dist=="WCCUSD" & rates$pre=="post"] -
              rates$rate_corr[rates$dist=="WCCUSD" & rates$pre=="pre"]))
  
  # corrected RD in control (post-pre)
  RD_cont = (rates$rate_corr[rates$dist=="WCCUSD" & rates$pre=="post"] -
               rates$rate_corr[rates$dist=="WCCUSD" & rates$pre=="pre"]) 
  
  # relative reduction adjusting for pre-int differencs
  did / RD_cont
  
}



##############################################
# Documentation:  correct_confounding
# Usage:          correct_confounding(corr_cases,
#                 RD_cd_pre, PCE1_pre, PCE0_pre,
#                 RD_cd_post, PCE1_post, PCE0_post, correct_hosp=T)
# Description:    calculate relative reduction difference-in-difference from 
#                 corrected cases counts 
#
# corr_cases:       data frame with corrected case counts 
# RD_cd_pre_prior:  prior distribution for risk difference for confounder & disease (pre-intervention)
# PCE1_pre_prior:   prior distribution for probability of confounder in intervention group (pre-intervention)
# PCE0_pre_prior:   prior distribution for probability of confounder in comparison group (pre-intervention)
# RD_cd_post_prior: prior distribution for risk difference for confounder & disease (post-intervention)
# PCE1_pre_prior:   prior distribution for probability of confounder in intervention group (post-intervention)
# PCE1_pre_prior:   prior distribution for probability of confounder in intervention group (post-intervention)
# correct_hosp:     boolean (T/F) to indicate whether to correct for flu testing rates among hospitalized  

# Returns: corrected relative reduction account for pre-intervention differences
# Output:  none
##############################################
# unmeasured confounding  ---------------------------------------
# a: count of disease in intervention group (E1)
# b: count of disease in comparison group (E0) 
# m: total in intervention group (E1)
# n: total in comparison group (E0)
# RD_cd: E(Y|C=1,A=0) - E(Y|C=0,A=0), assuming no EM 
# PCE1: P(C1|E1)
# PCE0: P(C1|E0)
# pre suffix is for pre-intervention period
# post suffix is for during intervention period
# formula based on Lash 2009 page 71
correct_confounding <- function(corr_cases,
                                RD_cd_pre, PCE1_pre, PCE0_pre,
                                RD_cd_post, PCE1_post, PCE0_post, correct_hosp=T){
  
  # use case counts corrected for incomplete flu testing
  if(correct_hosp){
    a_pre  = corr_cases$cases_corr[corr_cases$dist=="OUSD" & corr_cases$pre=="pre"]
    b_pre  = corr_cases$cases_corr[corr_cases$dist=="WCCUSD" & corr_cases$pre=="pre"]
    
    a_post  = corr_cases$cases_corr[corr_cases$dist=="OUSD" & corr_cases$pre=="post"]
    b_post  = corr_cases$cases_corr[corr_cases$dist=="WCCUSD" & corr_cases$pre=="post"]
    
    # use observed case counts
  }else{
    a_pre  = corr_cases$cases[corr_cases$dist=="OUSD" & corr_cases$pre=="pre"]
    b_pre  = corr_cases$cases[corr_cases$dist=="WCCUSD" & corr_cases$pre=="pre"]
    
    a_post  = corr_cases$cases[corr_cases$dist=="OUSD" & corr_cases$pre=="post"]
    b_post  = corr_cases$cases[corr_cases$dist=="WCCUSD" & corr_cases$pre=="post"]
  }
  
  m_pre  = corr_cases$N[corr_cases$dist=="OUSD" & corr_cases$pre=="pre"]
  n_pre  = corr_cases$N[corr_cases$dist=="WCCUSD" & corr_cases$pre=="pre"]
  
  m_post  = corr_cases$N[corr_cases$dist=="OUSD" & corr_cases$pre=="post"]
  n_post  = corr_cases$N[corr_cases$dist=="WCCUSD" & corr_cases$pre=="post"]
  
  m_pre  = corr_cases$N[corr_cases$dist=="OUSD" & corr_cases$pre=="pre"]
  n_pre  = corr_cases$N[corr_cases$dist=="WCCUSD" & corr_cases$pre=="pre"]
  
  m_post  = corr_cases$N[corr_cases$dist=="OUSD" & corr_cases$pre=="post"]
  n_post  = corr_cases$N[corr_cases$dist=="WCCUSD" & corr_cases$pre=="post"]
  
  a_pre_crude  = corr_cases$cases[corr_cases$dist=="OUSD" & corr_cases$pre=="pre"]
  b_pre_crude  = corr_cases$cases[corr_cases$dist=="WCCUSD" & corr_cases$pre=="pre"]
  
  a_post_crude  = corr_cases$cases[corr_cases$dist=="OUSD" & corr_cases$pre=="post"]
  b_post_crude  = corr_cases$cases[corr_cases$dist=="WCCUSD" & corr_cases$pre=="post"]
  
  # ---- pre intervention ----
  M1_pre = PCE1_pre * m_pre
  N1_pre = PCE0_pre * n_pre
  M0_pre = m_pre - M1_pre
  N0_pre = n_pre - N1_pre
  
  A1_pre = (M1_pre*(m_pre - M1_pre)*RD_cd_pre + M1_pre*a_pre)/m_pre
  B1_pre = (N1_pre*(n_pre - N1_pre)*RD_cd_pre + N1_pre*b_pre)/n_pre
  A0_pre = a_pre - A1_pre
  B0_pre = b_pre - B1_pre
  
  # ---- during intervention ----
  M1_post = PCE1_post * m_post
  N1_post = PCE0_post * n_post
  M0_post = m_post - M1_post
  N0_post = n_post - N1_post
  
  A1_post = (M1_post*(m_post - M1_post)*RD_cd_post + M1_post*a_post)/m_post
  B1_post = (N1_post*(n_post - N1_post)*RD_cd_post + N1_post*b_post)/n_post
  A0_post = a_post - A1_post
  B0_post = b_post - B1_post
  
  # ---- crude relative risks ----
  RD_crude_pre = (a_pre / m_pre) - (b_pre / n_pre)
  RD_crude_post = (a_post / m_post) - (b_post / n_post)
  
  # ---- crude relative reductions ----
  red_crude_pre = ((a_pre / m_pre) / (b_pre / n_pre)) - 1
  red_crude_post = ((a_post / m_post) / (b_post / n_post)) - 1
  
  # ---- DID ----
  DID_crude = ((a_post_crude / m_post) - (a_pre_crude / m_pre)) - 
    ((b_post_crude / n_post) - (b_pre_crude / n_pre))
  
  percent_red_DID_crude = (((a_post_crude / m_post) - (a_pre_crude / m_pre)) - 
                             ((b_post_crude / n_post) - (b_pre_crude / n_pre))) /  ((b_post_crude / n_post) - (b_pre_crude / n_pre))
  
  # pre-intervention mantel-haenszel risk difference
  num_pre = (((A1_pre * N1_pre - B1_pre * M1_pre)/(M1_pre + N1_pre))+
               ((A0_pre * N0_pre - B0_pre * M0_pre)/((M0_pre + N0_pre)))) 
  denom_pre =((N1_pre * M1_pre)/(N1_pre + M1_pre)) + 
    ((N0_pre * M0_pre)/(N0_pre + M0_pre))
  RD_mh_pre = num_pre / denom_pre
  RD_confound_pre = RD_crude_pre - RD_mh_pre
  
  # post-intervention mantel-haenszel risk difference
  num_post = (((A1_post * N1_post - B1_post * M1_post)/(M1_post + N1_post))+
                ((A0_post * N0_post - B0_post * M0_post)/((M0_post + N0_post)))) 
  denom_post =((N1_post * M1_post)/(N1_post + M1_post)) + 
    ((N0_post * M0_post)/(N0_post + M0_post))
  RD_mh_post = num_post / denom_post
  RD_confound_post = RD_crude_post - RD_mh_post
  
  # control group mantel-haenszel risk difference
  num_cont = ((B1_post*N1_pre - B1_pre*N1_post)/(N1_post + N1_pre)) +
    ((B0_post*N0_pre - B0_pre*N0_post)/(N0_post + N0_pre))
  denom_cont = (N1_post * N1_pre / (N1_post + N1_pre)) +
    (N0_post * N0_pre / (N0_post + N0_pre))
  
  RD_mh_cont = num_cont/denom_cont
  
  # bias-corrected DID
  DID_corrected = RD_mh_post - RD_mh_pre
  
  # DID due to confounding = DID crude - DID corrected
  DID_confounding = DID_corrected - DID_crude
  
  # relative reduction accounting for pre-intervention
  # differences and correcting for bias
  percent_red_DID_corrected = DID_corrected/RD_mh_cont
  
  RR_DID_confounding = (percent_red_DID_corrected+1)/(percent_red_DID_crude+1)
  red_crude_pre / percent_red_DID_corrected
  
  out = data.frame(DID_crude = DID_crude, 
                   DID_corrected = DID_corrected,
                   DID_confounding = DID_confounding,
                   percent_red_DID_crude = percent_red_DID_crude,
                   percent_red_DID_corrected = percent_red_DID_corrected,
                   RR_DID_confounding = RR_DID_confounding,
                   
                   RD_crude_pre = RD_crude_pre,
                   RD_crude_post = RD_crude_post,
                   RD_mh_pre = RD_mh_pre,
                   RD_mh_post = RD_mh_post,
                   RD_confound_pre = RD_confound_pre,
                   RD_confound_post = RD_confound_post,
                   
                   RD_cd_pre = RD_cd_pre, 
                   PCE1_pre = PCE1_pre, 
                   PCE0_pre = PCE0_pre,
                   RD_cd_post = RD_cd_post, 
                   PCE1_post = PCE1_post, 
                   PCE0_post = PCE0_post)
  
  return(out)
}

