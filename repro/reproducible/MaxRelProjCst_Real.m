%%
% MaxRelProjCst_Real.m
% Computes (a lower bound for) the real maximal relative projection constant
%
% Determines the quantity \lambda_R(m,N) described in the article
% "On Maximal Relative Projection Constants"
% by S. Foucart and L. Skrzypek
% based on the method presented in the supplementary note
%
% Usage: lambda = MaxRelProjCst_Real(m,N,nTest)
%
% m: positive integer representing the dimension of subspaces
% N: positive integer representing the dimension of the superspace
% nTest: positive integer representing the number of random initializations
%
% lambda: value of (a lower bound for) the real maximal relative projection constant

% Written by Simon Foucart in August 2016
% Send comments to simon.foucart@centraliens.net

function lambda = MaxRelProjCst_Real(m,N,nTest)

for t = 1:nTest
    [U_aux,~,~] = svd(randn(N,N));
    U = U_aux(:,1:m);
    [tt,~] = eigs(abs(U*U'),1,'LA');
    T = diag(tt);
    A = sign(U*U')';
    gamma_temp = 0;
    delta_temp = 1;
    while delta_temp - gamma_temp > 1e-7
        [U,evals] = eigs(T*A*T,m,'LA');
        gamma_temp = sum(sum(evals));
        gamma(t) = gamma_temp;
        [tt,delta_temp] = eigs(abs(U*U'),1,'LA');
        T = diag(tt);
        delta(t) = delta_temp;
        A = sign(U*U')';
    end
end

lambda = max(delta);

end