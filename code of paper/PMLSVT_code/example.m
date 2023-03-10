%% clear 

close all
clear all

%% loading file

Hours = xlsread('data.xlsx');
[m n] = size(Hours);

%% transform data based on weekdays to form a matrix with property of approximately low rank
day = 0;
M = zeros(24,105);
for i = 2 : m
    if ( Hours(i-1,2) ~= Hours(i,2) ) && (Hours(i,2)==6)
        day = day + 1;
    end
    M(Hours(i,1)+1,day) = Hours(i,3);
end

%% Establish the Bernoulli sampling 0-1 matrix and the corresponding observation matrix
% People can also use their own sampling 0-1 matrix by changing 'Omega'
[d1, d2] = size(M);
prop = 0.5;
Omega=binornd(1,prop,d1,d2);
Y = Omega.*M;


%% PMLSVT

%Setting parameters
K = 2000; lambda = 100; gamma = 1.1; L0 = 1e-04; beta = 1; alpha = 799;

%Call the main function
M1 = PMLSVT_completion(Y,Omega,alpha,beta,lambda, K,L0,gamma);

%% Display matrices
imagesc(Y); figure; imagesc(M1);