# Network meta-analysis of time-to-event data with Cox regression using individual patient-level data

This project contains the R functions for simulation and data analysis for NMA-IPD using stratified Cox model with random effects.

simulation.R
```
This code provides functions for generating NMA-IPD survival data under Weibull parametric proportional hazard models and comparisons from 4 NMA-IPD Cox model specifications:
1) (unstr_f) unstratified baseline hazard with fixed effect
2) (str_f) stratified baseline hazard with fixed effect
3) (unstr_r) unstratified baseline hazard with random effects
4) (str_r) stratified baseline hazard with random effects

The function "my.meta.survial.HR.sim()" requires the following arguments:
nt: number of clinical trial
ntrial: number of participants in each trial; vector of length nt
alphas: including mean log-HR comparing treatment k > 2 to the reference treatment (k = 1); vector of length K-1 (K is the number of treatments)
tau: variablity of the random effects; vector of length (K-1)
gamma: shape parameter of the baseline Weibull distribution; vector of length nt
lambda: scale parameter of the baseline Weibull distribution; vector of length nt
tcen: administrative censoring time
t_trt: treatments investigated in each trial; list of length nt

Values:

```
