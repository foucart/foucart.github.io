%%
% QMaxRelProjCst.m
% Computes the quasimaximal relative projection constant
%
% Determines the quantity \mu(m,N) described in the article
% "On Maximal Relative Projection Constants"
% by S. Foucart and L. Skrzypek
% based on the method presented in Section 3
%
% Usage: [mu,B,U] = QMaxRelProjCst(m,N)
%
% m: positive integer representing the dimension of subspaces
% N: positive integer representing the dimension of the superspace
%
% mu: the value of the quasimaximal relative projection constant
% B: NxN Seidel matrix yielding the maximum
% U: Nxm matrix yielding the maximum in the other expression for \mu

% Written by Simon Foucart in June 2015, modified in February 2016
% Send comments to simon.foucart@centraliens.net

function [mu,B,U] = QMaxRelProjCst(m,N)

% fetch all equivalent Seidel matrices, which have been precomputed
load('precomputed.mat','Seidel');
S = Seidel{N};    % each section of S is a Seidel matrix
K = size(S,3);    % number of equivalent Seidel matrices

% determine the maximum involving all the Seidel matrices,
% as well as the optimal matrix U
mu = 0;
I = eye(N);
for k = 1:K
    A = I + S(:,:,k);
    [evecs, evals] = eigs(A,m,'LA');
    mu_aux = sum(diag(evals))/N;
    if mu_aux > mu
        mu = mu_aux;
        B = A - I;
        U = evecs;
    end
end

end