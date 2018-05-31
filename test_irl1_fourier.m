%% IRL1 for Fourier incoherent design
close all
clear
clc

% we need CVX installed
try
    run('cvx_setup');
catch err
    error('CVX problem.');
end

% number of vectors
N = 29;
% dimension of the space
m = 8;

F = dftmtx(N)'; F = F(2:ceil(N/2)+1,:);

% how many trials
K = 50;
coherences = zeros(K, 1);
newsupports = zeros(m, K);
%% find K solutions
for k = 1:K
    if (k~=1)
        lambda = 10;
        w = ones(N, 1);
        epsi = 10e-7;
        zerosupport = randsample(2:N, ceil((N-m)/10));

        for i = 1:10
            cvx_begin
                cvx_quiet(true);
                variable x2(N)

                minimize norm(m^(-1)*F*x2, inf) + lambda/m*norm(w.*x2, 1)
                subject to
                         sum(x2) == m;
                         x2 <= 1;
                         x2 >= 0;
                         x2(zerosupport) == 0;
                         x2(1) == 1;
            cvx_end

            w = 1./(abs(x2) + epsi);
        end

        totalsupport = find(abs(x2)>10e-8);
        support = find(abs(x2)>10e-8);
    else
        support = (1:N)';
    end
    xnew = zeros(N, 1); xnew(support) = 1;

    % if there are more elements in the support, remove them
    for i = length(support):-1:m+1
        mincoherence = inf;
        for j = 1:length(support)
            newsupport = setdiff(support, support(j));
            thetry = zeros(N, 1); thetry(newsupport) = 1;
            newcoherence = norm(F*thetry, inf);

            if (newcoherence < mincoherence)
                mincoherence = newcoherence;
                index_found = j;
            end
        end
        support = setdiff(support, support(index_found));
        xnew = zeros(N, 1); xnew(support) = 1;
    end
    
    coherences(k) = norm(m^(-1)*F*xnew, inf);
    newsupports(:, k) = support;
end

%% accross all solutions, pick the best
mu = sqrt((N - m) / (m * (N - 1)));
[coherence, theindex] = min(coherences);
newsupport = newsupports(:, theindex);

save(['wl1 ' num2str(N) ' - ' num2str(m) '.mat'], 'newsupport', 'coherence');
