%% X-SIDCO TESTS

%% clean-up
close all
clear
clc

%% check if CVX is installed
try
    run('cvx_setup');
catch err
    error('CVX problem.');
end

%% experiment setup
for i = 1:1
    %% dimension of the space
    for n = 8
        %% number of vectors
        for N = 29
            %% number of iterations
            K = 500;

            %% performance limits
            umin1 = 1/sqrt(n);
            umin2 = sqrt((N-n)/(n*(N-1)));

            %% initialization
            A = randn(n, N);
            A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));

            %% improve initialization by SVD
            % real
            [U, ~, V] = svd(A);
            A = U*[eye(n) zeros(n, N-n)]*V';
            A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
            mcA = max(max(abs(A'*A) - eye(N)));
            
            % complex
            Ac = randn(n, N) + 1i*randn(n, N);
            Ac = bsxfun(@rdivide, Ac, sqrt(sum(abs(Ac).^2)));
            [U, ~, V] = svd(Ac);
            Ac = U*[eye(n) zeros(n, N-n)]*V';
            Ac = bsxfun(@rdivide, Ac, sqrt(sum(abs(Ac).^2)));
            mcAc = max(max(abs(Ac'*Ac) - eye(N)));
            
            % unital
            Au = Ac./abs(Ac);
            Au = bsxfun(@rdivide, Au, sqrt(sum(abs(Au).^2)));
            
            % real positive
            Ap = abs(A);
            
            % complex positive
            Acp = abs(randn(n, N)) + 1i*abs(randn(n, N));
            Acp = bsxfun(@rdivide, Acp, sqrt(sum(abs(Acp).^2)));

            %% X-SIDCO calls
            [Br, mcBr, timer, mcsBr] = r_sidco(A, K, umin2);
            [Bc, mcBc, timec, mcsBc] = c_sidco(Ac, K, umin2);

            [Brp, mcBrp, timerp, mcsBrp] = rp_sidco(Ap, K, umin2);
            [Bcp, mcBcp, timecp, mcsBcp] = cp_sidco(Acp, K, umin2);
            
            [Bu, mcBu, timeu, mcsBu] = u_sidco(Au, K, umin2);
            
            lambda = 1.8;
            [Bsr, mcBsr, timesr, mcsBsr] = sr_sidco(Br, K, lambda, umin2);
            [Bsc, mcBsc, timesc, mcsBsc] = sc_sidco(Bc, K, lambda, umin2);
            
            save(['incoherent frames ' num2str(N) ' - ' num2str(n) '.mat']);
        end
    end
end
