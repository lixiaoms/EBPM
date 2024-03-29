library(parallel)
source('mme_PLN.R')

source('update_A.R')

source('update_eta.R')

source('obj_gr_A.R')

source('est_A.R')

source('MAP.R')

source('est_A_obs.R')
source('obj_gr_A_obs.R')

## Main function for matrix denoising
EB_mle_denoise <- function(Y,k_max=100){
    a <- estA(Y)
    X_hat <- denoise_map2(a$mu,a$A,Y,obs=matrix(1,dim(Y)[1],dim(Y)[2]),k_max=k_max)
    return(X_hat)
}

## Main function for matrix completion
EB_mle_complete <- function(Y,obs,iter=20,k_max=100,control=0){
    Y <- Y*obs
    a <- estA_obs(Y,obs)
    X_hat <- denoise_map2(a$mu,a$A,Y,obs,k_max=k_max)
    if (iter<=0){
        ## Avoid the maximum value of estimated X is too big
        if (control==0){ return(X_hat) }
        return(pmin(X_hat,max(log(Y+0.5)*obs+1)))
    }
    NuNorm <- sum(diag(svd(X_hat)$d))
    for (i in 1:iter){
        ## Update Y by using observation data Y*obs and missing data entries of X_hat
        Y <- Y*obs+exp(X_hat)*(1-obs)
        a <- estA(Y)
        X <- denoise_map2(a$mu,a$A,Y,obs,k_max=k_max)
        ## Break if the nuclear norm of estimated X does not decrease
        NuNorm_1 <- sum(diag(svd(X)$d))
        if (NuNorm_1 > NuNorm*(1-1e-3)){
            return(pmin(X_hat,max(log(Y+0.5)*obs+0.5)))
        }
        ## Update X_hat
        X_hat <- X
        NuNorm <- NuNorm_1
    }
    ## Avoid the maximum value of estimated X is too big
    if (control==0){ return(X_hat) }
    return(pmin(X_hat,max(log(Y+0.5)*obs+0.5)))
}