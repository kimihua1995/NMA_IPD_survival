#a) genarating NMA-IPD survival data
my.meta.survial.HR.sim <- function(nt,betas,tau,gamma,lambda,alpha,t_trt){
  # nt: number of clinical trial
  # betas: including mean log-HR comparing treatment k > 2 to the reference treatment (k = 1); vector of length K-1 (K is the number of treatments arms)
  # tau: variablity of the random effects; vector of length (K-1)
  # gamma: shape parameter of the baseline Weibull distribution; vector of length nt
  # lambda: scale parameter of the baseline Weibull distribution; vector of length nt
  # alpha: values of trial-specific intercept; vector of length nt
  # t_trt: treatments investigated in each trial; list of length nt
  
   ntrial <- round(runif(nt,500,1000))
  ntrt <- length(betas) + 1
  betaj <- matrix(NA, nrow = nt, ncol = ntrt)
  betaj[,1] <- 0
  for (k in 2:ntrt){
    betaj[,k] <- rnorm(nt, mean = betas[k-1], sd = tau[k-1])
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
    timee <- (-log(u) * exp(-trt %*% betaj[j,] - alpha[j])/lambda[j])^(1/gamma[j])
    #timec <- runif(ntrial[j],0.1,tcen)
    timec <- pmin(rexp(ntrial[j],0.1),4)
    time <- pmin(timee,timec)
    #time[time > tcen] <- tcen
    #status <- as.numeric(timee < timec & timee < tcen)
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
dat <- my.meta.survial.HR.sim(nt = 4, 
                              betas = c(-0.2,-0.4), 
                              tau = c(0.1,0.2),
                              gamma = c(0.7,0.7,1.0,1.0), 
                              lambda = c(1,0.5,1,0.5),
                              alpha = rep(0,4),
                              t_trt = list(1:3,1:2,2:3,c(1,3)))
head(dat)






#b) Comparison of 9 Models in the simulation study
model_comp <- function(data){
  data$trt<-as.factor(1+data$trt2+2*data$trt3)
  data$trtbytrial<-as.factor(10*as.numeric(data$trial)+as.numeric(data$trt))
  # Model 1
  fit1 <- coxph(Surv(time,status) ~ trt2 + trt3 + trial,
                data = data)
  # Model 2
  fit2 <- coxme(Surv(time,status) ~ trt2 + trt3 + trial +
                  (trt2|trial) + (trt3|trial), data = data)
  # Model 3
  fit3 <- coxme(Surv(time,status) ~ trt2 + trt3 +
                  (1|trial), data = data)
  # Model 4
  fit4 <- coxme(Surv(time,status) ~ trt2 + trt3 + (1|trial) +
                  (trt2|trial) + (trt3|trial), data = data)
  # Model 5
  fit5 <- coxph(Surv(time,status) ~ trt2 + trt3 + strata(trial), 
                data = data)
  # Model 6
  fit6 <- coxme(Surv(time,status) ~ trt2 + trt3 + 
                  (trt2|trial) + (trt3|trial) + strata(trial), data = data)
  
  #############################
  # Models 7 to 9 use an alternative R argument to specify random treatment effects in the
  # coxme() function for our models. However, the between-study variability of the treatment 
  # effects are assumed to be equivalent under this argument.
  # The simulation results show that two arguments provide similar results. 
  # The results from Models 7 to 9 are not shown in the main manuscript.
  #############################
  # Model 7
  fit7 <- coxme(Surv(time,status) ~ trt + trial + (1|trtbytrial), data = data)
  # Model 8
  fit8 <- coxme(Surv(time,status) ~ trt + (1|trial/trt), data = data)
  # Model 9
  fit9 <- coxme(Surv(time,status) ~ trt + (1|trtbytrial) + strata(trial), data = data)
  
  coef.list <- list()
  coef.list[[1]] <- coef(fit1)[1:2]
  coef.list[[2]] <- coef(fit2)[1:2]
  coef.list[[3]] <- coef(fit3)[1:2]
  coef.list[[4]] <- coef(fit4)[1:2]
  coef.list[[5]] <- coef(fit5)[1:2]
  coef.list[[6]] <- coef(fit6)[1:2]
  coef.list[[7]] <- coef(fit7)[1:2]
  coef.list[[8]] <- coef(fit8)[1:2]
  coef.list[[9]] <- coef(fit9)[1:2]
  sd.list <- list()
  sd.list[[1]] <- sqrt(diag(vcov(fit1)))[1:2]
  sd.list[[2]] <- sqrt(diag(vcov(fit2)))[1:2]
  sd.list[[3]] <- sqrt(diag(vcov(fit3)))[1:2]
  sd.list[[4]] <- sqrt(diag(vcov(fit4)))[1:2]
  sd.list[[5]] <- sqrt(diag(vcov(fit5)))[1:2]
  sd.list[[6]] <- sqrt(diag(vcov(fit6)))[1:2]
  sd.list[[7]] <- sqrt(diag(vcov(fit7)))[1:2]
  sd.list[[8]] <- sqrt(diag(vcov(fit8)))[1:2]
  sd.list[[9]] <- sqrt(diag(vcov(fit9)))[1:2]
  
  
  AIC <- c(AIC(fit1),AIC(fit2),AIC(fit3),AIC(fit4),AIC(fit5),
           AIC(fit6),AIC(fit7),AIC(fit8),AIC(fit9))
  CI_L <- list()
  CI_H <- list()
  for (i in 1:9){
    CI_L[[i]] <- coef.list[[i]] - qnorm(0.975)*sd.list[[i]]
    CI_H[[i]] <- coef.list[[i]] + qnorm(0.975)*sd.list[[i]]
  }
  
  
  tau.list <- list()
  tau.list[[1]] <- NA
  tau.list[[2]] <- sqrt(unlist(fit2$vcoef))
  tau.list[[3]] <- NA
  tau.list[[4]] <- sqrt(unlist(fit4$vcoef))
  tau.list[[5]] <- NA
  tau.list[[6]] <- sqrt(unlist(fit6$vcoef))
  tau.list[[7]] <- sqrt(unlist(fit7$vcoef))
  tau.list[[8]] <- sqrt(unlist(fit8$vcoef))
  tau.list[[9]] <- sqrt(unlist(fit9$vcoef))
  
  return(list(coef=coef.list, sd=sd.list, 
              AIC=AIC,CI_L=CI_L, CI_H=CI_H, tau=tau.list))
}


#c) Simulation
my.sim <- function(seed,S,nt,betas,tau,t_trt,strata){
  # if strata == TRUE, set the simulation setting as stratified baseline hazards, meanwhile alpha = 0
  # if strata == FALSE, set the simulation setting as random trial-specific intercept, meanwhile alpha ~ N(0,1)

  
  coef.list <- rep(list(NULL),9)
  sd.list <- rep(list(NULL),9)
  tau.list <- rep(list(NULL),9)
  AIC.list <- NULL
  CI_L.list <- CI_H.list <- rep(list(NULL),9)
  
  set.seed(seed)
  for (i in 1:S){
    if (strata == TRUE){
      gamma <- runif(nt,0.5,1)
      lambda <- runif(nt,0.5,1)
      alpha <- rep(0,nt)
    }else if (strata == FALSE){
      gamma=rep(1,nt)
      lambda=rep(0.7,nt)
      alpha=rnorm(nt)
    }
    
    data <- my.meta.survial.HR.sim(nt = nt, betas = betas, tau = tau,
                                   gamma = gamma, lambda = lambda, 
                                   alpha = alpha, t_trt = t_trt)
    res <- model_comp(data)
    for (m in 1:9){
      coef.list[[m]] <- rbind(coef.list[[m]], res$coef[[m]])
      sd.list[[m]] <- rbind(sd.list[[m]], res$sd[[m]])
      CI_L.list[[m]] <- rbind(CI_L.list[[m]],res$CI_L[[m]])
      CI_H.list[[m]] <- rbind(CI_H.list[[m]],res$CI_H[[m]])
      tau.list[[m]] <- rbind(tau.list[[m]],res$tau[[m]])
    }
    AIC.list <- rbind(AIC.list, res$AIC)
    #if (i%%100 == 0) print(i)
  }
  
  beta_hat <- do.call("rbind",lapply(coef.list, colMeans))
  empri_sd <- do.call("rbind",lapply(coef.list, FUN = function(x){apply(x,2,sd)}))
  beta_sd <- do.call("rbind",lapply(sd.list, colMeans))
  #tau_hat <- do.call("rbind",lapply(tau.list, colMeans))
  CR <- NULL
  for (i in 1:9){
    CI_L <- t(t(CI_L.list[[i]]) <= betas)
    CI_H <- t(t(CI_H.list[[i]]) >= betas)
    CR <- rbind(CR, colMeans(CI_L * CI_H))
  }
  CR <- paste0(format(CR*100, nsmall = 1),"%")
  AIC.rank <- function(x) paste0(format(sum(x==1)/S*100,nsmall = 1),"%")
  AIC <- data.frame(AIC.mean = colMeans(AIC.list),
                    AIC.sd = apply(AIC.list,2,sd),
                    AIC.rank = apply(t(apply(AIC.list,1,rank)),2,AIC.rank))
  res = list(beta_hat = as.data.frame(beta_hat),
             bias = as.data.frame(beta_hat - rep(betas,each=9)),
             beta_sd = as.data.frame(beta_sd),
             empri_sd = as.data.frame(empri_sd),
             #tau_hat = as.data.frame(tau_hat),
             CR = data.frame(trt2=CR[1:9],trt3=CR[10:18]),
             AIC = AIC)

  return(list(res, tau.list))
}





t_trt=c(rep(list(1:3),5),rep(list(1:2),5),rep(list(c(1,3)),5))
betas=c(-0.2,-0.4)
nt=15

# Scenario 1
res_S1 <- my.sim(seed=12321, S=500, nt=nt, betas = betas, tau=c(0,0), 
                 t_trt = t_trt, strata = FALSE)

# Scenario 2
res_S2 <- my.sim(seed=12321, S=500, nt=nt, betas = betas, tau=c(0,0), 
                 t_trt = t_trt, strata = TRUE)

# Scenario 3
res_S3 <- my.sim(seed=12321, S=500, nt=nt, betas = betas, tau=c(0.1,0.2), 
                 t_trt = t_trt, strata = FALSE)

# Scenario 4
res_S4 <- my.sim(seed=12321, S=500, nt=nt, betas = betas, tau=c(0.1,0.2), 
                 t_trt = t_trt, strata = TRUE)
