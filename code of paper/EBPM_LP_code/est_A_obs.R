estA_obs <- function(Xdata,obs, mu.ini=NULL, A.ini=NULL,scale=NULL,max.iter=100){
    
       
  ##1. Preparation
  ##-------------------------------------------------------------------------------
  if (min(colSums(obs))==0){
      print("There is column with no observation, please delete it")
      break
      }
  if (min(rowSums(obs))==0){
      print("There is row with no observation, please delete it")
      break
  }
    
    Xdata <- as.matrix(Xdata*obs)
    n=nrow(Xdata)
    p=ncol(Xdata)

   if(is.null(mu.ini) || is.null(A.ini)){
        mme.value<- mmePLN(Xdata,obs)
        mu.array <- mme.value$mu
        a.array <- (1/mme.value$sigma)
        mu.old <- mu.array
        A.old <- diag(a.array)
    }
    if(!is.null(mu.ini)){mu.old = mu.ini}
    if(!is.null(A.ini)){A.old = A.ini}

    eig.tmp <- eigen(A.old, only.values=TRUE)
    alpha <- min(abs(eig.tmp$values),1)^2

    Amat.old <- A.old
    diff <- 1
    iter<- 0

    tmp <- objGrA_obs(Amat.old,Xdata,mu.old, obs, scale=scale)
    obj.old <- tmp$obj.ftn
    gr.old <- tmp$gr
    mu.old <- tmp$mu

    while(diff > 1e-6 && iter < max.iter){

        obj.diff<- -1
        if(alpha < 1e-3){
            eig.tmp <- eigen(Amat.old, only.values=TRUE)
            alpha <- min(abs(eig.tmp$values),1)^2
        }

        i.iter <- 0
        while(obj.diff < 0 && i.iter < 5){
            Amat.update <- updateA(Amat.old,gr.old,alpha) ## Checking if A.update <- positive-definite
            Amat.new <- Amat.update$A
            alpha <- Amat.update$alpha

            tmp <- objGrA_obs(Amat.new,Xdata,mu.old,obs, scale=scale)
            obj.new <- tmp$obj.ftn
            gr.new <- tmp$gr
            mu.new <- tmp$mu

            if(obj.new == -Inf){obj.new=-1e+100}
       #     if(obj.new==Inf){
      #          Amat.new <- Amat.old
       #         gr.new <- gr.old
      #          obj.new <- obj.old
       #         mu.new <- mu.old
        #    }

            obj.diff <- obj.old - obj.new

            if(is.na(obj.diff)){obj.diff <- 1}
            if(obj.diff <= 0){
                if(iter == 0){
                    alpha<- alpha/2
                }else{
                    alpha<- alpha/2
                    Amat.new <- Amat.old
                    gr.new <- gr.old
                    obj.new <- obj.old
                    mu.new <- mu.old
                }
            }
            i.iter <- i.iter +1
        }

        diff1 <- max(abs(Amat.old-Amat.new))/mean(abs(Amat.old))
        diff2 <- sum((Amat.old-Amat.new)^2)
        diff <- min(diff1,diff2)

        iter<- iter+1

        Amat.old <- Amat.new
        gr.old <- gr.new
        obj.old <- obj.new
        mu.old <- mu.new
    }

    tmp <- list(A = Amat.new, mu = mu.new, alpha= alpha, iter=iter)
    return(tmp)
}