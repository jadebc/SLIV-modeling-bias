# School-located influenza vaccination and community-wide indirect effects: comparing transmission models to epidemiologic models

## Overview
Mathematical models and empirical epidemiologic studies (e.g., randomized and observational studies) are complementary tools but may produce conflicting results for a given research question. We used sensitivity analyses and bias analyses to explore such discrepancies in a study of the indirect effects of influenza vaccination.

## Directory Structure

**`bias-analysis` :** folder containing probabilistic bias analysis scripts
* `0-base-functions` : R script containing general functions used across the analysis
* `0-config.R` : configuration file that sets data directories, sources base functions, and loads required libraries
* `1-did_percent_red_CI.R` : confidence interval for the percent reduction in influenza hospitalization 
* `2-define-hosp-priors.R` : define priors for bias analysis that corrects flu hospitalization case counts for incomplete influenza testing
* `3-hosp-cases.R` : correct hospitalization counts for incomplete flu testing same priors in each site
* `4-hosp-cases-diff-by-site.R` : correct hospitalization counts for incomplete flu testing allowing priors to differ between sites 
* `5-confounding.R` : correct DID and relative reduction for time-dependent confounding

**`mathematical-model` :** folder containing mathematical modeling scripts

**`figures` :** folder containing figure files
