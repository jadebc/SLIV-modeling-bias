###################################################
# Shoo the Flu evaluation bias analysis 

# confidence interval for the percent
# reduction in influenza hospitalization 

###################################################

rm(list=ls())

# load libraries, source file paths
source(paste0(here::here(), "/bias-analysis/0-config.R"))

# calculation by hand -------------------------------------------
d_full <- readRDS(data_path)
x=d_full$fluseasCDPH_2_5$eld

nt1 = sum(x$flucases[x$dist=="OUSD" & x$seas==1718])
nc1 = sum(x$flucases[x$dist=="WCCUSD" & x$seas==1718])
Nt1 = sum(x$N[x$dist=="OUSD" & x$seas==1718])
Nc1 = sum(x$N[x$dist=="WCCUSD" & x$seas==1718])

rate_t1 = nt1 / Nt1
rate_c1 = nc1 / Nc1

(rate_t1 - rate_c1)  / rate_c1

# DID functions -------------------------------------------

##############################################

# Documentation: fit_did_glm_seas
# Usage: fit_did_glm_seas(data, season, outcome, treatment, offset, covariates, label)
# Description: Preprocesses data and then fits a Poisson log-linear model 
#              and obtains an adjusted mean difference, DID and 95% CI
#              for a program year subtracting rates in pre-program years 
# Args/Options: 
# data:                    the data; must include a variable called seas, 
#                          indicating the program year / flu season
# season:                  number indicating program year / flu seas, as numeric
# preseas:                 pre-intervention seasons, as a numeric vector
# outcome:                 the outcome column name, as a string
# treatment:               the treatment column name, as a string   
# offset:                  log population offset column name, as a string
# covariates:              adjustment covariates, as a character vector
# label:                   description of dataset, as a string
# parameter:               "DID" for difference-in-differences; 
#                          "relDID" for relative reduction in DID
# 
# Returns: a formatted data frame that includes: adjusted difference estimate, 
# difference-in-difference estimate, standard error, the lower and upper bounds 
# of a 95% CI, and the season and person-years

fit_did_glm_seas = function(data, season, preseas, outcome, treatment, offset,
                            covariates, label, parameter = "DID"){
  pre_seas_list <- preseas
  analysis_seas_list <- c(pre_seas_list, season)
  
  fitdata_did <- data %>% 
    filter(seas %in% analysis_seas_list)
  
  fitdata_did$time = ifelse(fitdata_did$seas == season, 1, 0)
  
  cov_formula = paste(covariates, collapse = " + ")
  LHS = paste(outcome, "~")
  RHS_did = paste0(treatment, "*", "time")
  
  did_formula <- as.formula(paste(LHS, RHS_did))
  
  glm_did=glm(did_formula,offset=fitdata_did[[offset]],data=fitdata_did,
              family=poisson(link=log))
  
  res_relred = get_did_relred(glm_did)
  res_relred = res_relred %>% mutate(parameter = "Relative reduction difference-in-difference")
  
  N_did = fitdata_did %>% summarise(pyears = sum(N))
  
  res_relred = res_relred %>% mutate(pyears = N_did$pyears) 
  
  res = bind_rows(res_relred) %>% 
    mutate(seas = season,
           label = label,
           covariates = paste(covariates, collapse=", ")) %>%
    select(label, covariates, seas, pyears, everything()) %>%
    arrange(label, covariates, seas, parameter)
  
  return(res)
  
}




##############################################

# Documentation: get_did_relred
# Usage: get_did_relred(fit)
# Description: Fits a Poisson log-linear model and obtains the relative reduction DID and 95% CI
#              for a program year subtracting rates in pre-program years 
#              (rate_t1 - rate_c1) / rate_c1
# Args/Options: 
# fit: glm fit from a log-linear model with coefficients for treatment, pre/post
#      intervention, and an interaction term between the two and no other variables
# 
# Returns: a formatted data frame that includes: difference-in-difference estimate, 
# standard error, the lower and upper bounds of a 95% CI, and the season

get_did_relred = function(fit){
  
  b0 = fit$coefficients[1]
  b1 = fit$coefficients[2]
  b2 = fit$coefficients[3]
  b3 = fit$coefficients[4]
  
  reldid = (exp(b0+b1+b2+b3) - exp(b0+b2)) / exp(b0+b2)
  
  library(msm)
  g = as.formula(~ (exp(x1+x2+x3+x4) - exp(x1+x3)) / exp(x1+x3))
  se = deltamethod(g = g, mean = coef(fit), cov = vcov(fit))
  
  # 95% CI 
  lb = reldid - (qnorm(0.975)*se)
  ub = reldid + (qnorm(0.975)*se)
  
  # 2-sided p-value 
  t_stat = (reldid - 1) / se
  pval = 2 * pt(-abs(t_stat), df = nrow(fit$data) - 1)
  
  out = data.frame(estimate = reldid, se = se , lb = lb, ub = ub, pval = pval)
  rownames(out) = NULL
  
  return(out)
  
}




# DID relative reduction estimate and CI ---------------------------------------------
fit_did_glm_seas(
  data = x,
  season = 1718,
  preseas =  matrix(c(1112, 1213, 1314)),
  outcome = "flucases",
  treatment = "tr",
  offset = "logN",
  covariates = c("agecat", "sex", "race"),
  label = "Census data subset by zipcode (Primary)",
  parameter = "relDID"
)





