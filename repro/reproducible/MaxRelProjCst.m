%%
% MaxRelProjCst.m
% Computes the maximal relative projection constant
%
% Determines the quantity \lambda(m,N) described in the article
% "On Maximal Relative Projection Constants"
% by S. Foucart and L. Skrzypek
% based on the method presented in Section 3
%
% Usage: [lambda,B,t,U,V] = MaxRelProjCst(m,N,nIter)
%
% m: positive integer representing the dimension of subspaces
% N: positive integer representing the dimension of the superspace
% nIter: the number of iterations to be performed
%    (optional; a stopping criterion with a default tolerance is used otherwise)
%
% lambda: the value of (a lower bound for) the maximal relative projection constant
% B: NxN Seidel matrix yielding the maximum
% t: the unit vector yielding the maximum in either one of the expressions for \lambda
% U: Nxm matrix yielding the maximum in the alternative expression
% V: Nxm matrix whose columns span the subspace with maximal projection constant

% Written by Simon Foucart in June 2015, modified in February 2016
% Send comments to simon.foucart@centraliens.net

function [lambda,B,t,U,V] = MaxRelProjCst(m,N,nIter)

% fetch all equivalent Seidel matrices, which have been precomputed
load('precomputed.mat','Seidel');
S = Seidel{N};    % each section of S is a Seidel matrix
K = size(S,3);    % number of equivalent Seidel matrices
I = eye(N);
t = ones(N,1)/sqrt(N);    % initial choice of t^{(0)}

% determine the maximum involving all the Seidel matrices
% case 1: nIter is explicitely given
if nargin == 3
    for n = 1:nIter
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
    end
    
    % return the outputs not already assigned
    lambda = beta;       % value of the maximal projection constant
    V = diag(1./t)*U;    % subspace with maximal projection constant
end


% determine the maximum involving all the Seidel matrices
% case 2: nIter is not given
if nargin < 3
    TOL = 1e-6;        % arbitrary tolerance used as a stopping criterion
    NITERMAX = 50;     % arbitrary maximum number of iteration
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
    %  return the outputs not already assigned
    lambda = beta;       % value of the maximal projection constant
    V = diag(1./t)*U;    % subspace with maximal projection constant  
end

end