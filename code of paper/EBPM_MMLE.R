library(CVXR)
library(parallel)
library(PLNet)
library(orthopolynom)
library(poilog)

## Use mle_newton function in "PLNet" package to get sigma_12 in MLE of 2-dimensional Poisson-Log normal data
mmle_cov <- function(Y,obs,n=1,k_max=2){
    
    ## Parallel computing
    clus <- makeCluster(n)
    clusterExport(clus, "Y", envir=environment())
    clusterExport(clus, "obs", envir=environment())
    clusterExport(clus, "k_max", envir=environment())
    clusterExport(clus, "mle_newton", envir=environment())
    
    mmle_dim2 <- function(i){
        library(orthopolynom)
        output <- rep(0,ncol(Y))
        for (j in (i+1):ncol(Y)){
            index <- obs[,i]*obs[,j]==1
            if (sum(index)>7 && max(Y[index,i])>5 && max(Y[index,j])>5){
                cov_input_mle<-try(mle_newton(data_use = Y[index,c(i,j)],
                      S_depth = rep(1,sum(index)),
                      k_max = k_max,
                      core_num = 1),silent=TRUE)
                if('try-error' %in% class(cov_input_mle)){output[j]=0}
                else {output[j]<-cov_input_mle$mlesigmahat[[k_max+1]][1,2]}
                }
            }
        return(output)
    }
        
    cov_hat <- parApply(clus,matrix(1:(ncol(Y)-1),1,ncol(Y)-1),2,mmle_dim2)
    stopCluster(clus)
    
    cov_hat <- cbind(cov_hat,rep(0,ncol(Y)))
    cov_hat <- cov_hat+t(cov_hat)
    return(cov_hat)
    }

## Make estimated sigma positive semi-definite
sigma_semipositive<-function(cov_input){
  ## Find columns (rows) with all 0 and keep them unchange
  index <- apply(cov_input,1,function(x) all(x==0))
              
  ## Optimatization by minimizing the finite norm--------------------
  cov <- cov_input[!index,!index]
  S <- Variable(dim(cov)[1], dim(cov)[1], PSD = TRUE)
  
  obj <- Minimize(max(abs(S-cov)))
    
  prob <- CVXR::Problem(obj)
  
  result <- CVXR::solve(prob,solver="SCS",verbose=FALSE)
  sigma_me2 <- result$getValue(S)
                 
  cov_input[!index,!index] <- sigma_me2
  return(cov_input)
}

## Estimate X by maximizing posterior probability
## n is the core of parallel computing, k_max is the maximum iterations number in Newton's method
denoise_map <- function(mu,cov,Y,obs,n=1,k_max=10){
    Y[obs==0]=0
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
            X <- X-det_X
            if (max(abs(det_X))<1e-4){
                break
            }
        }
        return(X)
    }
    X_hat <- parApply(clus,data,1,column_newton)
    stopCluster(clus)
    return(t(X_hat))
    }

## Main function for matrix denoising
EB_mmle_denoise <- function(Y,n=1,k_max=10){
    obs <- matrix(1,dim(Y)[1],dim(Y)[2])
    ## Use mmle_cov function in "PLNet" package to get MMLE estimator of sigma and mu
    t=10
    cov_input_mle <- mle_newton(data_use = Y,
                      S_depth = rep(1,nrow(Y)),
                      k_max = t,
                      core_num = n)
    
    cov_mle <- sigma_semipositive(cov_input_mle$mlesigmahat[[t+1]])
    X_hat <- denoise_map(cov_input_mle$mlemu[[t+1]],cov_mle,Y,obs,n,k_max)
    return(X_hat)
}
                 
# Main function for matrix completion
EB_mmle_complete <- function(Y,obs,iter=3,n=1,k_max=10,control=0){
    Y <- Y*obs
    cov <- diag(ncol(Y))
    mu <- rep(0,ncol(Y))
    ## Use poilogMLE function in "poilog" package to get MLE of 1-dimensional Poisson-Log normal data (mu and sigma_jj)
    for (i in 1:ncol(Y)){
    log1 <- try(poilogMLE(Y[obs[,i]==1,i], startVals = c(mu=mu[i], sig=cov[i,i]),
    nboot = 0, zTrunc = 'False',
    method = "BFGS", control = list(maxit=100)),silent=TRUE)
    if('try-error' %in% class(log1)){log1 <- poilogMLE(Y[obs[,i]==1,i], startVals = c(mu=mu[i], sig=cov[i,i]),
    nboot = 0, zTrunc = 'False',
    method = "Nelder-Mead", control = list(maxit=100))}
    mu[i] <- log1$par[[1]]
    cov[i,i] <- log1$par[[2]]^2
    }
    
    ## Use mmle_cov function to get sigma_jk
    cov <- cov+mmle_cov(Y*obs,obs,n=n)
    
    cov_mle <- sigma_semipositive(cov)
    X_hat <- denoise_map(mu,cov_mle,Y,obs,n,k_max)
    if (iter<=0){
        ## Avoid the maximum value of estimated X is too big
        if (control==0){ return(X_hat) }
        return(pmin(X_hat,max(log(Y+0.5)*obs+1)))
    }
    for (i in 1:iter){
        ## Update Y by using observation data Y*obs and missing data entries of X_hat
        Y <- Y*obs+round(exp(X_hat))*(1-obs)
        t=10
        cov_input_mle <- mle_newton(data_use = Y,
                      S_depth = rep(1,nrow(Y)),
                      k_max = t,
                      core_num = n)
    
        cov_mle <- sigma_semipositive(cov_input_mle$mlesigmahat[[t+1]])
        X <- denoise_map(cov_input_mle$mlemu[[t+1]],cov_mle,Y,obs,n,k_max)
        
        ## Break if the maximum value of estimated X is too big
        if (max(X)>max(log(Y+0.5)*obs+0.5)){
            ## Avoid the maximum value of estimated X is too big
            if (control==0){ return(X_hat) }
            return(pmin(X_hat,max(log(Y+0.5)*obs+1)))
        }
        ## Update X_hat
        X_hat <- X
    }
    ## Avoid the maximum value of estimated X is too big
    if (control==0){ return(X_hat) }
    return(pmin(X_hat,max(log(Y+0.5)*obs+1)))
}