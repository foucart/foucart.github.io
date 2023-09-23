%%
% MaxRelProjCst_randomized.m
% This is basically the same code as MaxRelProjCst.m,
% the main difference being that the starting point t^{(0)}=[1,...,1]/sqrt(N)
% is changed to several random nonegative vectors of norm equal to 1.
% The number of random tests becomes an input,
% while the number of iterations is not a possible input anymore.
% The output is reduced to the maximum obtained over all random tests,
% and possibly the vector of all values obtained for all random tests.

% Written by Simon Foucart in February 2016
% Send comments to simon.foucart@centraliens.net

function [lambda,Lambda] = MaxRelProjCst_randomized(m,N,nTests)

% fetch all equivalent Seidel matrices, which have been precomputed
load('precomputed.mat','Seidel');
S = Seidel{N};    % each section of S is a Seidel matrix
K = size(S,3);    % number of equivalent Seidel matrices
I = eye(N);
TOL = 1e-6;        % arbitrary tolerance used as a stopping criterion
NITERMAX = 50;     % arbitrary maximum number of iteration

for test = 1:nTests
    t = rand(N,1); t = t/norm(t);    % initial choice of t^{(0)}
    % determine the maximum involving all the Seidel matrices
    n = 0;
    alpha = 0;
    beta = sqrt(N);
    while ( (n<NITERMAX) && (beta-alpha)>TOL )
        T = diag(t);
        % determination of \alpha_n and U_{(n)}
        alpha = 0;
        for k = 1:K
            A = I + S(:,:,k);
            [evecs,evals] = eigs(T*A*T,m,'LA');
            alpha_aux = sum(diag(evals));
            if alpha_aux > alpha
                alpha = alpha_aux;
                B = A - I;
                U = evecs;
            end
        end
        % determination of \beta_n and of t^{(n)}
        [t_aux,beta] = eigs(abs(U*U'),1,'LA');
        t = abs(t_aux);
        % increment n
        n = n+1;
    end
    Lambda(test) = beta;       % value of the maximal projection constant 
end

% return the output
lambda = max(Lambda);

end