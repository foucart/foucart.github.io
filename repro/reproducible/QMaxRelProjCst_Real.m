%%
% QMaxRelProjCst_Real.m
% Computes (a lower bound for) the real quasimaximal relative projection constant
%
% Determines the quantity \mu_R(m,N) described in the article
% "On Maximal Relative Projection Constants"
% by S. Foucart and L. Skrzypek
% based on the method presented in the supplementary note
%
% Usage: mu = QMaxRelProjCst_Real(m,N,nTest)
%
% m: positive integer representing the dimension of subspaces
% N: positive integer representing the dimension of the superspace
% nTest: positive integer representing the number of random initializations
%
% mu: value of (a lower bound for) the real quasimaximal relative projection constant

% Written by Simon Foucart in August 2016
% Send comments to simon.foucart@centraliens.net


function mu = QMaxRelProjCst_Real(m,N,nTest)

for t = 1:nTest
    [U_aux,~,~] = svd(randn(N,N));
    U = U_aux(:,1:m);
    A = sign(U*U')';
    alpha_temp = 0;
    beta_temp = 1;
    while beta_temp - alpha_temp > 1e-7
        [U,evals] = eigs(A,m,'LA');
        alpha_temp = sum(sum(evals))/N;
        alpha(t) = alpha_temp;
        beta_temp = sum(sum(abs(U*U')))/N;
        beta(t) = beta_temp;
        A = sign(U*U')';
    end
end

mu = max(beta);

end