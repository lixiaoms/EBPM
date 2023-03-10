#' updateA, a sub-function of estA()
#'
#' This function updates A as a sub-step of function, estA() and
#' checks whether A is positive definite
#'
#' @param A a pxp matrix
#' @param gr the gradient
#' @param alpha a step size
#' @return a list of A and alpha
updateA<- function( A, gr, alpha){
    
    gr <- (gr+t(gr))/2
    
    A.next <- A - alpha * gr
    
    check.diag <- sum(diag(A.next)<1e-5)
    check.eigen <- sum(eigen(A.next,only.values=TRUE)$values<1e-5)
    
    while( check.diag >0 || check.eigen >0){
        alpha <- alpha/2
        
        A.next <- A - alpha * gr
        check.diag <- sum(diag(A.next)<1e-5)
        check.eigen <- sum(eigen(A.next,only.values=TRUE)$values<1e-5)
    }
    
    tmp <- list(A = A.next, alpha = alpha)
    return(tmp)
}