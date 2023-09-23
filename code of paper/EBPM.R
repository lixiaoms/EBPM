library(CVXR)
library(parallel)

## The fist estimation of mu and sigma based on moment
moment_int<-function(data_use,obs){
  ##1. Preparation
  ##-------------------------------------------------------------------------------
  if (min(colSums(obs))==0){
      print("There is column with no observation, please delete it")
      break
  }
    
  if (min(colSums(data_use*obs))==0){
      print("There is column in which observed data are all 0, please delete it")
      break
  }
  ##-------------------------------------------------------------------------------
  
  ##2. Initialized estimation of mu and sigma
  ##-------------------------------------------------------------------------------
  data_obs<-data_use*obs  
  log_vec1<-as.vector(log(pmax(colSums(data_obs),0.5)/colSums(obs)))
  sigmahat<-t(log((t(data_obs) %*% data_obs)/(t(obs) %*% obs)) - log_vec1) - log_vec1
  ##
  diag(sigmahat)<-(log(colSums(data_obs * (data_obs - 1))/colSums(obs)) - 2 * log_vec1)
  sigmahat[!is.finite(sigmahat)]<-0
  mu<-log_vec1 - diag(sigmahat)/2
    return(list(mu = mu,
              sigmahat = sigmahat))
  
}

## The subsequent estimation of mu and sigma used in "estimate-complete" method
moment_sub<-function(data_use,obs){
  data_obs<-data_use*obs  
  log_vec1<-as.vector(log(pmax(colSums(data_use),0.5)/nrow(data_use)))
  sigmahat<-t(log((t(data_use) %*% data_use)/nrow(data_use)) - log_vec1) - log_vec1
  ##
  diag(sigmahat)<-(log(colMeans(data_use * data_use - data_obs)) - 2 * log_vec1)
  sigmahat[!is.finite(sigmahat)]<-0
  mu<-log_vec1 - diag(sigmahat)/2
    return(list(mu = mu,
              sigmahat = sigmahat))
  
}

## Make estimated sigma positive semi-definite
sigma_semipositive<-function(cov_input){
  ## Find columns (rows) with all 0 and keep them unchange
  index <- apply(cov_input,1,function(x) all(x==0))
              
  ## Optimatization by minimizing the infinity norm--------------------
  cov <- cov_input[!index,!index]
  S <- Variable(dim(cov)[1], dim(cov)[1], PSD = TRUE)
  
  obj <- Minimize(max(abs(S-cov)))
    
  prob <- CVXR::Problem(obj)
  
  result <- CVXR::solve(prob,solver="SCS",verbose=FALSE)
  sigma_me2 <- result$getValue(S)
                 
  cov_input[!index,!index] <- sigma_me2
  return(cov_input)
}

## Estimate X by maximizing posterior probability density
## n is the core of parallel computing, k_max is the maximum iterations number in Newton's method
denoise_map <- function(mu,cov,Y,obs,n=1,k_max=100){
    Y[obs==0]=0
    min_ev <- min(eigen(cov)$values)
    cov <- cov+diag(max(-min_ev,0),length(mu))
    ## Add a diagnonal matrix to make cov positive definite               
    cov <- cov+diag(1e-5,length(mu))
                   
    ## Get precision matrix
    prec <- solve(cov)
    ## Parallel computing
    clus <- makeCluster(n)
    clusterExport(clus, "mu", envir=environment())
    clusterExport(clus, "prec", envir=environment())
    clusterExport(clus, "k_max", envir=environment())
    ## MLE by newton method
    data <- cbind(Y,obs)
    column_newton <- function(data){
        len <- length(data)/2
        ## Split data to get Y and obs
        Y <- data[1:len]
        obs <- data[(len+1):(2*len)]
        X <- log(Y+0.5)
        for (k in 1:k_max){
            gradient <- prec%*%(X-mu)+(-Y+exp(X))*obs
            Hessian <- prec+diag(c(exp(X)*obs))
            det_X <- solve(Hessian)%*%gradient
            ## Line search
            value <- t(X-mu)%*%prec%*%(X-mu)/2+sum((-Y*X+exp(X))*obs)
            diff <- 1
            m <- 0
            while (diff>0 && m<100){
                X_1 <- X-(0.5)**m*det_X
                value_1 <- t(X_1-mu)%*%prec%*%(X_1-mu)/2+sum((-Y*X_1+exp(X_1))*obs)
                diff <- value_1-(value-0.5**m*0.1*t(gradient)%*%det_X)
                m <- m+1
            }
            ## Stop rule
            if (max(abs(X_1-X))<1e-4){
                X <- X_1
                break
            }
            X <- X_1
        }
        return(X)
    }
    X_hat <- parApply(clus,data,1,column_newton)
    stopCluster(clus)
    return(t(X_hat))
    }

## Main function for matrix denoising
EB_denoise <- function(Y,n=1,k_max=100,control=0){
    obs <- matrix(1,dim(Y)[1],dim(Y)[2])
    moment <- moment_int(Y,obs)
    cov_moment <- sigma_semipositive(moment$sigmahat)
    X_hat <- denoise_map(moment$mu,cov_moment,Y,obs,n,k_max)
    ## whether to avoid the maximum value of estimated X is too big
    if (control==0){ return(X_hat) }
    return(pmin(X_hat,max(log(Y+0.5)*obs+1)))
}

## Main function for matrix completion
EB_complete <- function(Y,obs,iter=20,n=1,k_max=100,control=0){
    Y[is.na(Y)] <- 0
    moment <- moment_int(Y,obs)
    cov_moment <- sigma_semipositive(moment$sigmahat)
    X_hat <- denoise_map(moment$mu,cov_moment,Y,obs,n,k_max)
    if (iter<=0){
        ## whether to avoid the maximum value of estimated X is too big
        if (control==0){ return(X_hat) }
        return(pmin(X_hat,max(log(Y+0.5)*obs+1)))
    }
    NuNorm <- sum(diag(svd(X_hat)$d))
    for (i in 1:iter){
        ## Update Y by using observation data Y*obs and missing data entries of X_hat
        Y <- Y*obs+exp(X_hat)*(1-obs)
        moment <- moment_sub(Y,obs)
        cov_moment <- sigma_semipositive(moment$sigmahat)
        X <- denoise_map(moment$mu,cov_moment,Y,obs,n,k_max)
        ## Break if the nuclear norm of estimated X does not decrease
        NuNorm_1 <- sum(diag(svd(X)$d))
        if (NuNorm_1 > (NuNorm-1e-3)*(1-1e-3)){
            if (control==0){ return(X_hat) }
            return(pmin(X_hat,max(log(Y+0.5)*obs+1)))
        }
        ## Update X_hat
        X_hat <- X
        NuNorm <- NuNorm_1
    }
    ## whether to avoid the maximum value of estimated X is too big
    if (control==0){ return(X_hat) }
    return(pmin(X_hat,max(log(Y+0.5)*obs+1)))
}