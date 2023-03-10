#' etaUpdate
#'
#' This function computes eta.star for given mu and A to maximize the likelihood
#'
#' @param mu a mean of p genes, px1 vector
#' @param A an inverse of covariance matrix, pxp matrix
#' @param Xn counts for a sample, px1
#' @param eta.old initial estimate of eta
#' @param scale.n a K vector for a sample, px1 vector
#' @param eps a criterion for convergence where d(eta.old, eta.new) < eps
#' @param max.iter the maximum number of iterations
#' @return a list of eta.star and sigma.star
#' @export
etaUpdate <- function(mu, A, Xn, index=rep(TRUE,length(mu)), scale.n = NULL,eps= 1e-6, max.iter=100){
    
    if (sum(index)==1){
    mu <- mu[index]
    A <- A[index,index]
    Xn <- Xn[index]
    
    if(is.null(scale.n)){scale.n <- 1}
   eta.new <- log(Xn+1)
   diff <- 1
   iter<-1
   
   while(diff>eps && iter <= max.iter){
      eta <- eta.new
      grad <- Xn-(scale.n * exp(eta))-(eta-mu)*A
      tmp.inv <- (scale.n*exp(eta)+A)^(-1)
     
      eta.new <- eta + tmp.inv * grad
      diff<- abs(eta-eta.new)
      iter<-iter+1
   }
    
   tmp <- list(eta.star = eta.new, sigma.star = tmp.inv)
   return(tmp)}
    
    else{
    mu <- mu[index]
    A <- A[index,index]
    Xn <- Xn[index]
        
    if(is.null(scale.n)){scale.n <- (Xn*0 + 1)}
   eta.new <- log(Xn+1)
   p <- length(Xn)
   diff <- 1
   iter<-1
   
   while(diff>eps && iter <= max.iter){
      eta <- c(eta.new)
      grad <- t(Xn-(scale.n * exp(eta))-(eta-mu)%*%A)
      tmp.inv <- solve(diag(scale.n*c(exp(eta)))+A)
     
      eta.new <- eta + tmp.inv %*% grad
      diff<- max(abs(eta-eta.new))
      iter<-iter+1
   }
    
   tmp <- list(eta.star = eta.new, sigma.star = tmp.inv)
   return(tmp)
}
}