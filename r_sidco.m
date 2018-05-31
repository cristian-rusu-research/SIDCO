function [bestA, minmc, time, mcs] = r_sidco(A, K, umin)
%%% R-SIDCO
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
A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
mcs = max(max(abs(A'*A) - eye(m)));
bestA = A;
minmc = mcs;

%% the norm of the constraint
nrmCnst = 2;

%% main iterations
oldA = A;
ordering = randperm(m);

for k=1:K
    
    for i = ordering
        target = A(:, i);
        Awork = A(:, [1:i-1 i+1:end]);
        
        %% change signs to simplify optimization problem
        signs = (Awork'*target <= 0);
        for j=1:m-1
            if (signs(j) == 1)
                Awork(:, j) = -Awork(:,j);
            end
        end
        
        corre = Awork'*target;
        T = max(corre);
        
        %% main optimization problem
        cvx_begin
            variable Anew(n,1);

            cvx_quiet(true);
            % the max sometimes gives problems to the solver, especially
            % when m is close to n
            %minimize ( max(Awork'*Anew) )
            minimize ( norm(Awork'*Anew, inf) )
            subject to
                    norm(Anew-target, nrmCnst) <= sqrt(1-T^2);
        cvx_end

        %% normalize solution
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
            A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
            ordering = randperm(m);
        end
    end

    %% check if limit was reached
    if (abs(mc - umin) < 10e-5)
        break;
    end
end

time = toc;
