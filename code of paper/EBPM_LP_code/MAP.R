denoise_map2 <- function(mu,prec,Y,obs,n=1,k_max=10){
    Y[obs==0]=0
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
