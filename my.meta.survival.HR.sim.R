# genarating NMA-IPD survival data
my.meta.survial.HR.sim <- function(nt,ntrial,alphas,tau,gamma,lambda,tcen,t_trt){
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
dat <- my.meta.survial.HR.sim(nt = 4, 
                              ntrial = round(runif(4,500,1000)),
                              alphas = c(-0.2,-0.4), 
                              tau = c(0.01,0.01),
                              gamma = c(0.7,0.7,1.0,1.0), 
                              lambda = c(1,1.5,1,1.5),
                              tcen = 2, 
                              t_trt = list(1:3,1:3,1:3,1:3))
head(dat)
