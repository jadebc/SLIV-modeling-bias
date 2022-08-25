###################################################
# Shoo the Flu evaluation bias analysis 

# correct hospitalization counts 
# for incomplete flu testing
# allowing priors to differ between sites 

# not a probabilistic bias analysis

# output: corrected case counts
###################################################

rm(list =ls())

# load libraries, source file paths
source(paste0(here::here(), "/bias-analysis/0-config.R"))

# load data ---------------------------------------
d_full <- readRDS(data_path)
d <- d_full$fluseasCDPH_2_5

eld <- d$eld

all_priors <- read.csv(paste0(here::here(), "/results/priors_hosp.csv"))

# define priors ---------------------------------------
draw_prior_dis_65 = sample(all_priors$prior_dis_65, size = 1)

unif_prior = runif(n = 10000,
                   min = min(all_priors$prior_dis_65),
                   max = max(all_priors$prior_dis_65))

prior_list = seq(min(all_priors$prior_dis_65),
                      max(all_priors$prior_dis_65), by=0.01)

# create data frame with all combinations of sequence 
# of priors by each site
prior_df = expand.grid(prior_ousd = prior_list, 
                       prior_wcc = rev(prior_list))

# calculate corrected relative reduction ---------------------------------------
prior_df$relred = NA
for(i in 1:nrow(prior_df)){
  prior_df$relred [i] <- data.frame(relred = get_did_rates_dist(data = d$eld, 
                                     ousd_prior = prior_df$prior_ousd[i], 
                                     wcc_prior = prior_df$prior_wcc[i]) %>% unlist()) 
}
prior_df$relred=unlist(prior_df$relred)
prior_df$col =   ifelse(prior_df$relred <0 & prior_df$relred> -0.02, 
                        "#067BEC", "#C8C8C8")
prior_df$match =   as.factor(ifelse(prior_df$relred <0 & prior_df$relred> -0.02, 
                        "Match transmission model", "Does not match"))
prior_df = prior_df %>% 
  mutate(prior_diff = prior_ousd - prior_wcc)

# make a plot ---------------------------------------
panel_a <- ggplot(prior_df , aes(x = prior_ousd*100, y = prior_wcc*100)) + 
  geom_point(aes(col = match)) +
  scale_color_manual(values = c("#C8C8C8", "#37C17B")) + 
  xlab("% tested for influenza\n(intervention)") + 
  ylab("% tested for influenza\n(comparison)") + 
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  ggtitle("A)")

panel_b <- ggplot(prior_df, aes(x=prior_diff*100, y=relred*100)) + 
  geom_point(aes(col = prior_ousd*100), alpha=0.25) + 
  ylab("Relative reduction in\n influenza hospitalizations") + 
  xlab("Difference in % tested for\ninfluenza between sites") + 
  scale_color_continuous("% tested in\nintervention site", high = "#132B43", low = "#56B1F7")+
  theme_minimal() +
  theme(legend.position = "bottom")+
  ggtitle("B)")

plot <- grid.arrange(panel_a, panel_b, ncol=2)
ggsave(plot, filename = paste0(here::here(), "/figures/fig-flu-testing.png"),
                               width=7, height=4)

# summarize model match relative reductions ---------------------------------------
# proportion of results that match transmission model 
prop.table(table(prior_df$match))
min(prior_df$prior_ousd[prior_df$match=="Match transmission model"])
max(prior_df$prior_ousd[prior_df$match=="Match transmission model"])

min(prior_df$prior_diff[prior_df$relred>0])
max(prior_df$prior_diff[prior_df$relred>0])



