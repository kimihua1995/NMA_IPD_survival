#a) genarating NMA-IPD survival data
NMAIPD.survival.dat <- function(nt,ntrial,alphas,tau,gamma,lambda,tcen,t_trt){
  # nt: number of clinical trial
  # ntrial: number of participants in each trial; vector of length nt
  # alphas: including mean log-HR comparing treatment k > 2 to the reference treatment (k = 1); vector of length K-1 (K is the number of treatments arms)
  # tau: variablity of the random effects; vector of length (K-1)
  # gamma: shape parameter of the baseline Weibull distribution; vector of length nt
  # lambda: scale parameter of the baseline Weibull distribution; vector of length nt
  # tcen: administrative censoring time
  # t_trt: treatments investigated in each trial; list of length nt
  
  ntrt <- length(alphas) + 1
  alphaj <- matrix(NA, nrow = nt, ncol = ntrt)
  alphaj[,1] <- 0
  for (k in 2:ntrt){
    alphaj[,k] <- rnorm(nt, mean = alphas[k-1], sd = tau[k-1])
  }
  
  dat_comb <- list()
  
  for (j in 1:nt){
    n_trt_j <- length(t_trt[[j]])
    prob <- rep(0,ntrt)
    prob[t_trt[[j]]] <- 1
    trt0 <- sample(x = 1:ntrt, ntrial[j], replace = T, prob = prob)
    trt <- matrix(0,nrow = ntrial[j],ncol = ntrt)
    for (i in 1:ntrial[j]){
      trt[i,trt0[i]] <- 1
    }
    
    u <- runif(ntrial[j],0.00001,1)
    timee <- (-log(u) * exp(-trt %*% alphaj[j,])/lambda[j])^(1/gamma[j])
    timec <- runif(ntrial[j],0.1,tcen)
    time <- pmin(timee,timec)
    status <- as.numeric(timee < timec)
    dat <- as.data.frame(cbind(time,status,trt))
    names(dat) <- c("time","status",paste0("trt",1:ntrt))
    dat$trial <- j
    dat_comb[[j]] <- dat
  }
  
  data <- do.call("rbind", dat_comb)
  data$trial <- factor(data$trial)
  return(data)
}



# sample code
set.seed(12321)
dat <- NMAIPD.survival.dat(nt = 4, 
                           ntrial = round(runif(4,500,1000)),
                           alphas = c(-0.2,-0.4), 
                           tau = c(0.01,0.01),
                           gamma = c(0.7,0.7,1.0,1.0), 
                           lambda = c(1,1.5,1,1.5),
                           tcen = 2, 
                           t_trt = list(1:3,1:2,2:3,c(1,3)))
head(dat)






#b) model comparison
model_comp <- function(data){
  # data: output from NMAIPD.survival.dat()
  require(survival)
  require(coxme)
  fit1 <- coxph(Surv(time,status) ~ trt2 + trt3, data = data)
  fit2 <- coxph(Surv(time,status) ~ trt2 + trt3 + strata(trial), data = data)
  fit3 <- coxme(Surv(time,status) ~ trt2 + trt3 + 
                  (trt2|trial) + (trt3|trial), data = data)
  fit4 <- coxme(Surv(time,status) ~ trt2 + trt3 + 
                  (trt2|trial) + (trt3|trial) + strata(trial), data = data)
  
  coef.list <- list()
  coef.list[[1]] <- coef(fit1)
  coef.list[[2]] <- coef(fit2)
  coef.list[[3]] <- coef(fit3)
  coef.list[[4]] <- coef(fit4)
  sd.list <- list()
  sd.list[[1]] <- sqrt(diag(vcov(fit1)))
  sd.list[[2]] <- sqrt(diag(vcov(fit2)))
  sd.list[[3]] <- sqrt(diag(vcov(fit3)))
  sd.list[[4]] <- sqrt(diag(vcov(fit4)))
  tau.list <- list()
  tau.list[[1]] <- sqrt(unlist(fit3$vcoef))
  tau.list[[2]] <- sqrt(unlist(fit4$vcoef))
  AIC <- c(AIC(fit1),AIC(fit2),AIC(fit3),AIC(fit4))
  CI_L <- list()
  CI_H <- list()
  for (i in 1:4){
    CI_L[[i]] <- coef.list[[i]] - qnorm(0.975)*sd.list[[i]]
    CI_H[[i]] <- coef.list[[i]] + qnorm(0.975)*sd.list[[i]]
  }
  
  return(list(coef=coef.list, sd=sd.list, tau=tau.list, AIC=AIC,
              CI_L=CI_L, CI_H=CI_H))
}


model_comp(dat)


