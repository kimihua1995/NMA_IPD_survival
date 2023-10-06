# PH assumption checking and extended Cox models
# example code for covariate unadjusted models only
# can be naturally extended to covariate adjusted model
# we use simulated data from NMAIPD.survival.dat()
set.seed(12321)
data <- my.meta.survial.HR.sim(nt = 4, 
                              betas = c(-0.2,-0.4), 
                              tau = c(0.1,0.2),
                              gamma = c(0.7,0.7,1.0,1.0), 
                              lambda = c(1,0.5,1,0.5),
                              alpha = rep(0,4),
                              t_trt = list(1:3,1:2,2:3,c(1,3)))

library(coxme)
library(coxph)
# 1) PH assumption
# using Trial 2 and trt 2 vs. trt 1 as our example

## G1: log cumulaitve hazard plots
library(survminer)
legend.labs <- c("trt1","trt2")
palette <- c("#7570B3","#E7298A")
fit_G1 <- survfit(Surv(time,status) ~ trt2,data = data[data$trial==2,])
curve <- ggsurvplot(fit,data = data,size = 0.5,
                    censor.size = 1,conf.int = F,
                    pval = F,
                    xlim = c(0.001,2),
                    xlab = "Log time",
                    ylab = "Log cumulative hazard",
                    title = "Trial 2",
                    surv.plot.height = 5,
                    risk.table = F,
                    legend = "bottom",
                    legend.labs = legend.labs,
                    legend.title = "",
                    palette = palette,
                    fontsize = 3.0,
                    fun = "cloglog")
curve

## G2: scaled Schoenfeld's residual plots
fit_G2 <- coxph(Surv(time,status) ~ trt2, 
              data = data[data$trial==2,])
ylab <- "trt2"
main <- "trt2 vs. trt1 (Trial 2)"
plot(cox.zph(fit_G2), var = 1, resid = F, ylim = c(-2,2),
     xlab = "Time",ylab = ylab, lwd = 2,
     main = main)
abline(h=0, col = "red", lty = "dotted", lwd = 2)
abline(h=coef(fit_G2)[1], col = "darkblue", lty = "longdash", lwd = 2)
legend("topright", c("Null","Cox"),
       col = c("red","darkblue"),
       lty = c("dotted","longdash"),
       bty = 'n', border = NA)


## G3: plots of observed cumulative residual with simulated process under PH assumption
fit_G3 <- cox.aalen(Surv(time,status) ~ prop(trt2), 
            weighted.test = 0, data = data[data$trial==2,])
plot(fit_G3,score = T,main = c("trt2 vs. trt1 (Trial 2)"))



## T1: add treatment-by-lot(t) interactions
fit_T1 <- coxme(Surv(time,status) ~ trt2 +  trt2:log(time), data = data[data$trial==2,])
summary(fit_T1)


## T2: scaled Schoenfeld's residuals score test
cox.zph(fit_G2)

## T3: goodness-of-fit test on cumulative residuals
fit_G3








# 2) Extended Cox models
## a) piecewise Cox model
## t1, t2: two time cutpoints
data_split <- survSplit(Surv(time, status) ~ ., data = data,
                       cut = c(t1,t2), episode = "tgroup", id = "id")
data_split$p1 <- as.numeric(data_split$tgroup == 1)
data_split$p2 <- as.numeric(data_split$tgroup == 2)
data_split$p3 <- as.numeric(data_split$tgroup == 3)

fit_piece <- coxme(Surv(tstart,time,status) ~ trt2:p1 + trt2:p2+ trt2:p3+ 
               trt3:p1 + trt3:p2+ trt3:p3+
               (trt2|trial) + (trt3|trial) + strata(trial), data = data)


## b) cumulative Cox model
## t_cum: the selected cumulative time
data$status_cum <- data$status
data$time_cum <- data$time
data$status_cum[data$time > t_cum] <- 0
data$time_cum[data$time > t_cum] <- t_cum

fit_cum <- coxme(Surv(time_cum,status_cum) ~ trt2 + trt3 + 
          (trt2|trial) + (trt3|trial) + strata(trial), data = data)








