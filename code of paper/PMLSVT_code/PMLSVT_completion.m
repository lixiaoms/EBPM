
function [ M1 ] = PMLSVT_completion( Y, Omega, alpha, beta, lambda, K, L0, gamma)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Given the matrix with missing Poisson observations, the function will recovery 
%the original matrix by nuclear norm regularization optimization.   


%Y : d1-by-d2 observation matrix with missing entries. 
%Omega : d1-by-d2 0-1 matrix that denotes the location of observed entries.
%alpha : Upper bound for the estimator
%beta : Lower bound for the estimator
%lambda : Balance parameter. Large lambda will return smooth result and small lambda will return spike result
%K : Largest number of iterations 
%L0 : Initial reciprocal of step size.
%gamma : Multiplier for the update of L0

%More details about the parameters can be found in the following paper: 
%Cao, Yang, and Yao Xie. "Poisson Matrix Recovery and Completion." arXiv preprint arXiv:1504.05229 (2015).


%Author: Yang Cao
%Date: 11/12/2015
%Version: 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    %% initialization
    [d1, d2] = size(Y);
    M0 = Y;
    for i = 1 : d1
        for j = 1 : d2
            if (M0(i,j)==0)
                M0(i,j) = (alpha+beta)/2;
            end
        end
    end
    
    %% PMLSVT_completion
    for k = 1 : K
        p0 = sum(sum(Omega.*(M0-Y.*log(M0+1e-15))))+lambda*sum(svd(M0));
        Gradient = Omega.*(1 - Y./(M0+1e-15));
        while (1)
            C = M0 - 1/L0*Gradient;
            [U,D,V] = svd(C,'econ');
            D = diag(max(diag(D)-lambda/L0,0));
            M1 = U*D*V';
            M1 = max(beta,M1);
            M1 = min(alpha,M1);
            p1 = sum(sum(Omega.*(M1-Y.*log(M1+1e-15))))+lambda*sum(diag(D));
            minus = p1-p0;
            if ( (minus<0) || (abs(minus)<1e-05) )
                break;
            end
            L0 = L0 * gamma;
        end
        M0=M1;
    end


end

