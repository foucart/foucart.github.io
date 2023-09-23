%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reproducible MATLAB file accompanying the paper
% A LOWER BOUND FOR THE FIFTH MAXIMAL PROJECTION CONSTANT
% by B. Deregowska, M. Fickus, S. Foucart, and B. Lewandowska
% 
% Written in May 2022
% Send comments to simon.foucart@centraliens.net
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;


%% The upper and lower bounds for the maximal relative projection constant
% the upper bound, del, is always valid,
% the lower bound, gam, is only valid under the existence of certain MUETFs

del = @(m,n) m/n*(1+sqrt((n-1)*(n-m)/m));
gam = @(m,k,l) (m - del(m,k)*del(m,l))/(2*sqrt(m)-del(m,k)-del(m,l));


%% E.g., the lower bound evaluated at (5,6,10) is approximately 2.06919

fprintf(strcat('The lower bound for lambda(5,10) is', 32,...
    num2str(gam(5,6,10),'%.5f'),'\n'))


%% Verify that the construction of Section 3 produces MUETFs
% note: the function MUETF is defined at the end of the file

s = 3;
[V,W,m,k,l] = MUETF(s);

fig1=figure(1);
subplot(4,2,1)
imagesc(V*V')
colorbar
title('Representation of VV^T to show that V is a tight frame')
xlabel(strcat('k/m=',num2str(k/m)))
subplot(4,2,2)
imagesc(abs(V'*V))
colorbar
title('Representation of |V^TV| to show that V is equiangular')
xlabel(strcat('\phi_{m,k}=',num2str(sqrt((k-m)/m/(k-1)))))
subplot(4,2,3)
imagesc(W*W')
colorbar
title('Representation of WW^T to show that W is a tight frame')
xlabel(strcat('l/m=',num2str(l/m)))
subplot(4,2,4)
imagesc(abs(W'*W))
colorbar
title('Representation of |W^TW| to show that W is equiangular')
xlabel(strcat('\phi_{m,l}=',num2str(sqrt((l-m)/m/(l-1)))))
subplot(4,2,[5:8])
imagesc(abs(V'*W))
colorbar
title('Representation of |V^TW| to show that V and W are mutually unbiased')
xlabel(strcat('1/m^{1/2}=',num2str(1/sqrt(m))))
fig1.WindowState = 'maximized';


%% Verify that the candidates t=t_\theta and U=U_\theta associated with MUETFs 
% do satisfy the two necessary conditions for maximality

x_star = (del(m,k)-del(m,l))/(2*sqrt(m)-del(m,k)-del(m,l));
theta = acos(x_star)/2;
t = [cos(theta)/sqrt(k)*ones(k,1); sin(theta)/sqrt(l)*ones(l,1)];
U = [cos(theta)*sqrt(m/k)*V sin(theta)*sqrt(m/l)*W];
% first necessary condition
[ldg_evec,~] = eigs(abs(U'*U),1,'largestreal');
sgn = t(1)/ldg_evec(1);
fprintf(strcat(...
    'The leading eigenvector of |U^T*U| equals t (2-norm of difference=',...
    num2str(norm(t-sgn*ldg_evec)),')\n'))
% second necessary condition
T = diag(t);
evals = sort(eig(T*sign(U'*U)*T),'descend');
ldg_evals = evals(1:m);
fprintf(strcat(...
    'The rows of U are m leading eigenvectors of T*sgn(U^T*U)*T (Frobenius-norm of difference=',...
    num2str(norm(T*sign(U'*U)*T*U' - U'*diag(ldg_evals),'fro')),')\n'))


%% Interesting observation (proved in the article):
% the m leading eigenvalues of T*sgn(U'*U)*T are all equal

fprintf(strcat('The 1-st and m-th largest eigenvalues of T*sgn(U^T*U)*T are equal to',...
    32, num2str(ldg_evals(1)), ' and to ', 32, num2str(ldg_evals(m)),'\n'))


%% In the case m=21, k=28, l=36, there is a better lower bound for lambda(m)
% than gamma(21,28,36)---lambda(21,126) equals its upper bound delta(21,126) 

fprintf(strcat('For lambda(21), the lower bound gamma(21,28,36)=',...
    num2str(gam(21,28,36),'%.5f'), ' is not as good as the lower bound delta(21,126)=',...
    num2str(del(21,126),'%.5f'),'\n'))


%% In the case m=5, k=6, l=10, there does not seems to be a better lower bound
% for lambda(5) than gamma(5,6,10)

m = 5;
% experiment with several values of N, n_test, n_iter
% e.g N=20, n_test=500, n_iter=2000: approx 3 min,
% and N=32, n_test=500, n_iter=2000: approx 7 min.
N = 25;
n_test = 500;
n_iter = 2000;
% the following estimate lambda(m,N) from below by alternating maximization 
res = zeros(1,n_test);
for i=1:n_test
    [U,~,~] = svd(randn(N,m),'econ');
    U = U';
    [t,~] = eigs(abs(U'*U),1,'largestreal');
    T = diag(t);
    A = sign(U'*U);
    for j=1:n_iter
        [UU,~] = eig(T*A*T);
        U = UU(:,N-m+1:N)';
        [t,res_temp] = eigs(abs(U'*U),1,'largestreal');
        T = diag(t);
        A = sign(U'*U);   
    end
    res(i) = res_temp;
end
% comparison of the best alternating maximization lower bound
% with the gamma(5,6,10), conjectured to be the value of lambda(5)
fprintf(strcat('The best lower bound for lambda(5,',num2str(N),') is', 32,...
    num2str(max(res),'%.5f'),' <=', 32, num2str(gam(5,6,10),'%.5f'),'\n'))

%% the function realizing the construction of MUETFs from Section 3

function [V,W,m,k,l] = MUETF(s) 

m = (2^(2*s)-1)/3;              % dimension of the frame vectors
k = 2^(s-1)*(2^s-1);            % number of vectors in the first frame
l = 2^(s-1)*(2^s+1);            % number of vectors in the second frame
Y = dec2bin(0:2^(2*s)-1)-'0';   % all binary vectors of length 2s
Q = kron(eye(s),[0 1;0 0]);     % kernel of quadratic form on F_2^(2s)
B = kron(eye(s),[0 1;1 0]);     % kernel of bilinear form on F_2^(2s)
Gamma = (-1).^(mod(Y*B*Y',2));  % corresponding character table of F_2^(2s)
DS = mod(diag(Y*Q*Y'),2);
D  = find(DS==1);               % first of the paired difference sets 
DC = find(DS==0);               % second of the paired difference sets 
X = Gamma(D,DC);                % the matrix obeying conditions (i)-(iv)
[P,~,Q] = svds(X,m);            % the compact SVD of X
V = sqrt(k/m)*P';               % first of the mutually unbiased ETFs
W = sqrt(l/m)*Q';               % second of the mutually unbiased ETFs

end
