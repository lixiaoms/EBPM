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
