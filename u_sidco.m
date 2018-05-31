function [bestA, minmc, time, mcs] = u_sidco(A, K, umin)
%%% U-SIDCO
%% check if CVX is installed
try
    run('cvx_setup');
catch err
    error('CVX problem.');
end

%% start timer
tic;

%% setup
[n, m] = size(A);

%% initialization
A = A./abs(A);
A = bsxfun(@rdivide, A, sqrt(sum(abs(A).^2)));
mcs = max(max(abs(A'*A) - eye(m)));
bestA = A;
minmc = mcs;

%% the norm of the constraint
nrmCnst = 2;

%% main iterations
ordering = randperm(m); 
oldA = A;

for k=1:K
    
    for i = ordering
        target = A(:, i);
        Awork = A(:, [1:i-1 i+1:end]);
        
        %% main optimization problem
        cvx_precision medium;
        cvx_solver('sedumi');
        cvx_begin
            variable Anew(n,1) complex;

            cvx_quiet(true);
            minimize ( max(abs(Awork'*Anew)) )
            subject to
                    norms(Anew-target, nrmCnst, 2) <= 0.1*1/sqrt(n);
                    norms(Anew, 1, 2) - 1/sqrt(n) <= 0.05*1/sqrt(n);
        cvx_end

        %% normalize solution
        Anew = Anew./abs(Anew);
        Anew = Anew/norm(Anew);
        A(:, i) = Anew;
        oldA(:, i) = target;
    end
    
    %% new frame, new mutual coherence
    mc = max(max(abs(A'*A) - eye(m)));

    %% update best
    if (mc < minmc)
        bestA = A;
        minmc = mc;
    end
    
    %% keep track of coherences
    mcs = [mcs mc];
    
    %% check convergence
    if (k >= 5)
        if (std(mcs(end-2:end)) < 10e-6)
            [U, ~, V] = svd(A);
            A = U*[eye(n) zeros(n, m-n)]*V';
            A = A./abs(A);
            A = bsxfun(@rdivide, A, sqrt(sum(abs(A).^2)));
        end
    end
    
    %% change ordering
    ordering = randperm(m);
    
    %% check if limit was reached
    if (abs(mc - umin) < 10e-5)
        break;
    end
end

time = toc;
