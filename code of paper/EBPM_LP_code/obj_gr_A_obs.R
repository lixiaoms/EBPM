objGrA_obs <- function(A, Xdata, mu,obs, scale=NULL){
    
    if(is.null(scale)){ scale = Xdata*0+1 }
    
    p<-ncol(Xdata)
    n<-nrow(Xdata)
    
    
    gr.sum1 <- matrix(0,p,p)
    gr.sum2 <- 0
    log.l.sum <- 0
    
    ## mu is updated by newton's method
    mu.gr <- mu*0
    mu.he <- matrix(0,p,p)
    
    for(i in 1:n){
        index <- obs[i,]==1
        
        if (sum(index)==1){
            
        scale.n <- scale[i,index]
        eta.tmp <- etaUpdate(mu,A,Xdata[i,],index=index, scale.n=scale.n)
        eta.i.star <- eta.tmp$eta.star
        sigma.i.star <- eta.tmp$sigma.star
        
        mu.gr[index] <- mu.gr[index] + A[index,index] * (eta.i.star - mu[index])
        mu.he[index,index] <- mu.he[index,index] + A[index,index]
        
        eta.mu <- (eta.i.star-mu[index])
        gr.i <- -(A[index,index])^(-1) + eta.mu^2 + sigma.i.star
        gr.sum1[index,index] <- gr.sum1[index,index] + gr.i
        
        A.star <- scale.n * exp(eta.i.star) + A[index,index]
        avg.diag <- A.star
        A.star.bar <- A.star/avg.diag
            
        log.l.i <- 0.5*log(A[index,index]) + (-1)*0.5*(eta.i.star-mu[index])* A[index,index] *(eta.i.star-mu[index])- scale.n*exp(eta.i.star) + eta.i.star * Xdata[i,index] - 0.5*(log(A.star.bar)+ sum(index)*log(avg.diag))}
        
        else{
        scale.n <- scale[i,index]
        eta.tmp <- etaUpdate(mu,A,Xdata[i,],index=index, scale.n=scale.n)
        eta.i.star <- eta.tmp$eta.star
        sigma.i.star <- eta.tmp$sigma.star
            
        mu.gr[index] <- mu.gr[index] + A[index,index] %*% (eta.i.star - mu[index])
        mu.he[index,index] <- mu.he[index,index] + A[index,index]
        
        eta.mu <- (eta.i.star-mu[index])
        gr.i <- -solve(A[index,index]) + eta.mu %*% t(eta.mu) + sigma.i.star
        gr.sum1[index,index] <- gr.sum1[index,index] + gr.i
        
        A.star <- diag(scale.n * c(exp(eta.i.star))) + A[index,index]
        avg.diag <- mean(diag(A.star))
        A.star.bar <- A.star/avg.diag
            
        log.l.i <- 0.5*log(det(A[index,index])) + (-1)*0.5*t(eta.i.star-mu[index])%*% A[index,index] %*%(eta.i.star-mu[index])- sum(scale.n*exp(eta.i.star)) + sum(eta.i.star * Xdata[i,index])- 0.5*(log(det(A.star.bar))+ sum(index)*log(avg.diag))
            }
        log.l.sum <- log.l.sum + log.l.i
    }
    gr.mean <- gr.sum1/n
    gr.mean <- gr.mean/2
    log.l.mean <- log.l.sum/n
    mu.mean <- mu + solve(mu.he) %*% mu.gr
    
    gr.A <- gr.mean
    obj.A <- (-1)*log.l.mean
    
    tmp <- list(obj.ftn = obj.A, gr = gr.A, mu = c(mu.mean))
    return(tmp)
}
