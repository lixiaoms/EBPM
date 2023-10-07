%% PMLSVT
maxNumCompThreads(1);

%Setting parameters
K = 2000; lambda = 100; gamma = 1.1; L0 = 1e-04; beta = 0.1; alpha = 1000;


%% change the length of matrix
q=100; r=10;
for p = 100:100:1000
    t1 = []; t2 = [];
    for n = 1:10
        Y = readmatrix(strcat('..\simulated data\',string(p),'-100-10-',string(n),'.txt'));
        obs = readmatrix(strcat('..\simulated data\',string(p),'-100-10-',string(n),'-obs.txt'));
        
        tic
        M = PMLSVT_completion(Y,ones(p,q),alpha,beta,lambda, K,L0,gamma);
        t1 = [t1, toc];
        save(strcat('..\PMLSVT\',string(p),'-100-10-',string(n),'-denoise.txt'), 'M', '-ascii');
        
        tic
        M = PMLSVT_completion(Y.*obs,obs,alpha,beta,lambda, K,L0,gamma);
        t2 = [t2, toc];
        save(strcat('..\PMLSVT\',string(p),'-100-10-',string(n),'-complete.txt'), 'M', '-ascii');
    end
    save(strcat('..\PMLSVT\',string(p),'-100-10-denoise time.txt'), 't1', '-ascii');
    save(strcat('..\PMLSVT\',string(p),'-100-10-complete time.txt'), 't2', '-ascii');
end

%% change the rank of matrix
p=500; q=100;
for r = 5:5:50
    for n = 1:10
        Y = readmatrix(strcat('..\simulated data\500-100-',string(r),'-',string(n),'.txt'));
        obs = readmatrix(strcat('..\simulated data\500-100-',string(r),'-',string(n),'-obs.txt'));
        
        M = PMLSVT_completion(Y,ones(p,q),alpha,beta,lambda, K,L0,gamma);
        save(strcat('..\PMLSVT\500-100-',string(r),'-',string(n),'-denoise.txt'), 'M', '-ascii');
        
        M = PMLSVT_completion(Y.*obs,obs,alpha,beta,lambda, K,L0,gamma);
        save(strcat('..\PMLSVT\500-100-',string(r),'-',string(n),'-complete.txt'), 'M', '-ascii');
    end
end

%% change the proportion of observations
p=500; q=100;
for r = 0.1:0.1:0.9
    for n = 1:10
        Y = readmatrix(strcat('..\simulated data\500-100-10-',string(n),'.txt'));
        obs = readmatrix(strcat('..\simulated data\500-100-',string(n),'-',string(r),'-obs.txt'));
        
        M = PMLSVT_completion(Y.*obs,obs,alpha,beta,lambda, K,L0,gamma);
        save(strcat('..\PMLSVT\500-100-',string(n),'-',string(r),'-complete.txt'), 'M', '-ascii');
    end
end

%% HSI
% p=400; q=191;
% t = [];
% for n = 1:25
%     Y = readmatrix(strcat('..\simulated data\HSI\',string(n),'.txt'));
%     tic
%     M = PMLSVT_completion(Y,ones(p,q),alpha,beta,lambda, K,L0,gamma);
%     t = [t, toc];
%     save(strcat('..\PMLSVT\HSI\',string(n),'-denoise.txt'), 'M', '-ascii');
% end
% save(strcat('..\PMLSVT\HSI\denoise time.txt'), 't', '-ascii');
% 
% for r = 0.1:0.1:0.9
%     t=[];
%     for n = 1:25
%         Y = readmatrix(strcat('..\simulated data\HSI\',string(n),'.txt'));
%         obs = readmatrix(strcat('..\simulated data\HSI\',string(n),'-',string(r),'-obs.txt'));
%         tic
%         M = PMLSVT_completion(Y.*obs,obs,alpha,beta,lambda, K,L0,gamma);
%         t = [t, toc];
%         save(strcat('..\PMLSVT\HSI\',string(n),'-',string(r),'-complete.txt'), 'M', '-ascii');
%     end
%     save(strcat('..\PMLSVT\HSI\',string(r),'-complete time.txt'), 't', '-ascii');
% end