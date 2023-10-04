# Network meta-analysis with individual participant-level data of time-to-event outcomes using Cox regression

This project contains the R functions for simulation and data analysis for NMA-IPD

simulation.R
```
This code provides functions for generating NMA-IPD survival data under Weibull parametric proportional hazard models and comparisons from 4 NMA-IPD Cox model specifications:

a) function "my.meta.survial.HR.sim()"
  Arguments:
  # nt: number of clinical trial
  # betas: including mean log-HR comparing treatment k > 2 to the reference treatment (k = 1); vector of length K-1 (K is the number of treatments arms)
  # tau: variablity of the random effects; vector of length (K-1)
  # gamma: shape parameter of the baseline Weibull distribution; vector of length nt
  # lambda: scale parameter of the baseline Weibull distribution; vector of length nt
  # alpha: values of trial-specific intercept; vector of length nt
  # t_trt: treatments investigated in each trial; list of length nt


  Values: 
  a data frame including:
  time: time to event
  status: censoring indicator; 1 if event, 0 if censored
  trt1 -- trtK: treatment indicator, K is the number of treatment arms
  trial: trial number from 1 to nt




b) function "model_comp()"
9 IPD-NMA unadjusted Cox model specifications are compared:
Model 1: stratified Cox model (Equation 1) with fixed treatment effect.
Model 2: Cox model (Equation 2) using fixed trial-specific intercept and fixed treatment effect.
Model 3: Cox model (Equation 2) using random trial-specific intercept and fixed treatment effect.
Model 4: stratified Cox model (Equation 1) with random treatment effect.
Model 5: Cox model (Equation 2) using fixed trial-specific intercept and random treatment effect (uising specification ).
Model 6: Cox model (Equation 2) using random trial-specific intercept and random treatment effect.
Model 7: stratified Cox model (Equation 1) with random treatment effect (specification a).
Model 8: Cox model (Equation 2) using fixed trial-specific intercept and random treatment effect (specification a).
Model 9: Cox model (Equation 2) using random trial-specific intercept and random treatment effect (specification a).

require packages "survival" and "coxme"

Arguments:
data: output from from my.meta.survial.HR.sim()

Values:
a list including:
coef: point estimation of coefficients in the cox models
sd: standard deviation estimation of coefficients in the cox models
tau: estimation of between-trial variability
AIC: AIC of each model
CI_L, CI_H: lower and upper bound of 95% confidence interval for coef
```


analysis.R
```
This file contains the code used for PH assumption assessment and extended Cox models for data analysis.
We only use the covariate unadjusted model only as our example, which can be easily extended to the covarate adjusted models
The apprroaches for PH assumption assessment includes:
G1: log cumulaitve hazard plots
G2: scaled Schoenfeld's residual plots
G3: plots of observed cumulative residual with simulated process under PH assumption
T1: add treatment-by-lot(t) interactions
T2: scaled Schoenfeld's residuals score test
T3: goodness-of-fit test on cumulative residuals

The extended Cox models are:
a) piecewise Cox model
b) cumulative Cox model

```

