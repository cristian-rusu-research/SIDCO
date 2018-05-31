function [bestA, minmc, time, mcs, pz] = sr_sidco(A, K, lambda, umin)
%%% SR-SIDCO
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

for k=1:K
    
    for i = 1:m
        target = A(:, i);
        Awork = A(:, [1:i-1 i+1:end]);
        
        signs = (Awork'*target <= 0);
        for j=1:m-1
            if (signs(j) == 1)
                Awork(:, j) = -Awork(:,j);
            end
        end
        
        corre = Awork'*target;
        T = max(corre);
        
        %% main optimization problem
        cvx_precision medium;
        cvx_solver('sedumi');
        cvx_begin
            variable Anew(n,1);

            cvx_quiet(true);
            minimize ( max(Awork'*Anew) + 1/n*lambda*norm(Anew,1))
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
    if (k > 10)
        if (std(mcs(end-10:end)) < 10e-5)
            break;
        end
    end
    
    %% check if limit was reached
    if (abs(mc - umin) < 10e-5)
        break;
    end
end


%% second set of iterations
oldA = A;
minmc = inf;
for k=1:100
    
    tostop = 0;
    for i = 1:m
        target = A(:, i);
        Awork = A(:, [1:i-1 i+1:end]);
        
        signs = (Awork'*target <= 0);
        for j=1:m-1
            if (signs(j) == 1)
                Awork(:, j) = -Awork(:,j);
            end
        end
        
        corre = Awork'*target;
        notsupport = find(abs(target)<10e-5);
        T = max(corre);
        
        %% main optimization problem
        cvx_precision medium;
        cvx_solver('sedumi');
        cvx_begin
            variable Anew(n,1);

            cvx_quiet(true);
            minimize ( max(Awork'*Anew) )
            subject to
                    norm(Anew-target, nrmCnst) <= sqrt(1-T^2);
                    Anew(notsupport) == 0;
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
    if (k > 10)
        if (std(mcs(end-10:end)) < 10e-5)
            break;
        end
    end
    
    %% check if limit was reached
    if (abs(mc - umin) < 10e-5)
        break;
    end
end

pz = length(find(abs(bestA)<10e-5))/(n*m)*100;
time = toc;
