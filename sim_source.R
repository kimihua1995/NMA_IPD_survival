# data generation
my.meta.survial.HR.sim <- function(nt,betas,tau,gamma,lambda,tcen,t_trt){
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
    timee <- (-log(u) * exp(-trt %*% betaj[j,])/lambda[j])^(1/gamma[j])
    timec <- runif(ntrial[j],0.1,tcen)
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




# model compare
model_comp <- function(data){
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




# simulation function
my.sim <- function(seed,nt,betas,tau,gamma.list,lambda.list,tcen,t_trt){
  coef.list <- list(NULL,NULL,NULL,NULL)
  sd.list <- list(NULL,NULL,NULL,NULL)
  tau.list <- list(NULL,NULL)
  AIC.list <- NULL
  CI_L.list <- CI_H.list <- list(NULL,NULL,NULL,NULL)
  
  set.seed(seed)
  for (i in 1:10000){
    #gamma <- sample(gamma.list,nt,replace = T)
    #lambda <- sample(lambda.list,nt,replace = T)
    gamma <- runif(nt,gamma.list[1],gamma.list[2])
    lambda <- runif(nt,lambda.list[1],lambda.list[2])
    data <- my.meta.survial.HR.sim(nt = nt, betas = betas, tau = tau,
                                   gamma = gamma, lambda = lambda, tcen = tcen, t_trt = t_trt)
    res <- model_comp(data)
    for (m in 1:4){
      coef.list[[m]] <- rbind(coef.list[[m]], res$coef[[m]])
      sd.list[[m]] <- rbind(sd.list[[m]], res$sd[[m]])
      CI_L.list[[m]] <- rbind(CI_L.list[[m]],res$CI_L[[m]])
      CI_H.list[[m]] <- rbind(CI_H.list[[m]],res$CI_H[[m]])
      if (m < 3){
        tau.list[[m]] <- rbind(tau.list[[m]],res$tau[[m]])
      }
    }
    AIC.list <- rbind(AIC.list, res$AIC)
    if (i%%1000 == 0) print(i)
  }
  
  beta_hat <- do.call("rbind",lapply(coef.list, colMeans))
  empri_sd <- do.call("rbind",lapply(coef.list, FUN = function(x){apply(x,2,sd)}))
  beta_sd <- do.call("rbind",lapply(sd.list, colMeans))
  tau_hat <- do.call("rbind",lapply(tau.list, colMeans))
  CR <- NULL
  for (i in 1:4){
    CI_L <- t(t(CI_L.list[[i]]) <= betas)
    CI_H <- t(t(CI_H.list[[i]]) >= betas)
    CR <- rbind(CR, colMeans(CI_L * CI_H))
  }
  CR <- paste0(format(CR*100, nsmall = 1),"%")
  AIC.rank <- function(x) paste0(format(sum(x==1)/10000*100,nsmall = 1),"%")
  AIC <- data.frame(AIC.mean = colMeans(AIC.list),
                    AIC.sd = apply(AIC.list,2,sd),
                    AIC.rank = apply(t(apply(AIC.list,1,rank)),2,AIC.rank))
  res = list(beta_hat = as.data.frame(beta_hat),
             bias = as.data.frame(beta_hat - rep(betas,each=4)),
             beta_sd = as.data.frame(beta_sd),
             empri_sd = as.data.frame(empri_sd),
             tau_hat = as.data.frame(tau_hat),
             CR = data.frame(trt2=CR[1:4],trt3=CR[5:8]),
             AIC = AIC)

  return(res)
}
