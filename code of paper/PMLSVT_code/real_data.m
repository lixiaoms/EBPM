%% PMLSVT
%% maxNumCompThreads(1);

%% hic data
K = 2000; lambda = 100; gamma = 1.1; L0 = 1e-04; beta = 1; alpha = 100;
for r = 0.5:0.1:0.9
    for n = 1:10
        Y = readmatrix(strcat('..\real data\hic_chr22_24_32_36mb.csv'));
        obs = readmatrix(strcat('..\real data\hic-',string(n),'-',string(r),'-obs.txt'));
        M = PMLSVT_completion(Y.*obs,obs,alpha,beta,lambda, K,L0,gamma);
        save(strcat('..\PMLSVT\hic-',string(r),'-',string(n),'-complete.txt'), 'M', '-ascii');
    end
end

%% bike data
K = 2000; lambda = 100; gamma = 1.1; L0 = 1e-04; beta = 1; alpha = 799;
for r = 0.5:0.1:0.9
    for n = 1:10
        Y = readmatrix(strcat('..\real data\bike.csv'));
        obs = readmatrix(strcat('..\real data\bike-',string(n),'-',string(r),'-obs.txt'));
        M = PMLSVT_completion(Y.*obs,obs,alpha,beta,lambda, K,L0,gamma);
        save(strcat('..\PMLSVT\bike-',string(r),'-',string(n),'-complete.txt'), 'M', '-ascii');
    end
end

%% brain image data
K = 2000; lambda = 100; gamma = 1.1; L0 = 1e-04; beta = 1; alpha = 5000;
for r = 0.5:0.1:0.9
    for n = 1:10
        Y = readmatrix(strcat('..\real data\brain_image.csv'));
        obs = readmatrix(strcat('..\real data\brain_image-',string(n),'-',string(r),'-obs.txt'));
        M = PMLSVT_completion(Y.*obs,obs,alpha,beta,lambda, K,L0,gamma);
        save(strcat('..\PMLSVT\brain_image-',string(r),'-',string(n),'-complete.txt'), 'M', '-ascii');
    end
end