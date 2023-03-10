#' objGrA, sub-function of estA()
#'
#' This function computes the objective function and its gradient for given mu, A, Z, U, W and Xdata
#'
#' @param A a pxp concentration matrix
#' @param Xdata a datat nxp matrix where n is the number of samples and p is the number of genes
#' @param mu mean of p genes, px1 vector
#' @param scale a nxp matrix of scale K
#' @return a list of obj.A, gr.A and updated mu
#' @export
objGrA <- function(A, Xdata, mu, scale=NULL){
    
    if(is.null(scale)){ scale = Xdata*0+1 }
    
    p<-ncol(Xdata)
    n<-nrow(Xdata)
    
    inv.A <- solve(A)
    det.A <- det(A)
    diag.A <- diag(A)
    
    gr.sum1 <- (-1)*n*inv.A
    gr.sum2 <- 0
    log.l.sum <- n*0.5*log(det.A)
    mu.tmp <- mu*0
    
    for(i in 1:n){
        scale.n <- scale[i,]
        eta.tmp <- etaUpdate(mu,A,Xdata[i,],scale.n=scale.n)
        eta.i.star <- eta.tmp$eta.star
        sigma.i.star <- eta.tmp$sigma.star
        
        mu.tmp <- mu.tmp + t(eta.i.star)
        
        eta.mu <- (eta.i.star-mu)
        gr.i <- eta.mu %*% t(eta.mu) + sigma.i.star
        gr.sum1 <- gr.sum1 + gr.i
        
        A.star <- diag(scale.n * c(exp(eta.i.star))) + A
        avg.diag <- mean(diag(A.star))
        A.star.bar <- A.star/avg.diag
        
        log.l.i <- (-1)*0.5*t(eta.i.star-mu)%*% A %*%(eta.i.star-mu)- sum(scale.n*exp(eta.i.star)) + sum(eta.i.star * Xdata[i,])- 0.5*(log(det(A.star.bar))+ p*log(avg.diag))
        log.l.sum <- log.l.sum + log.l.i
    }
    gr.mean <- gr.sum1/n
    gr.mean <- gr.mean/2
    log.l.mean <- log.l.sum/n
    mu.mean <- mu.tmp/n
    
    gr.A <- gr.mean
    obj.A <- (-1)*log.l.mean
    
    tmp <- list(obj.ftn = obj.A, gr = gr.A, mu = c(mu.mean))
    return(tmp)
}
