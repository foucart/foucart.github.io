%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reproducible file containing the experiments mentioned in the article
% "On Maximal Relative Projection Constants"
% by S. Foucart and L. Skrzypek
% Created by Simon Foucart in June 2015
% Last modified in August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%% Computation of the Real QuasiMaximal Relative Projection Constants \mu(m,N)
% the function QMaxRelProjCst(m,N) is defined in its separate file

Nmax = 10;
TableMu = NaN(Nmax,Nmax);      % table to contain the values of \mu
TableMuB = cell(Nmax,Nmax);    % table to contain the optimal Seidel matrices
TableMuU = cell(Nmax,Nmax);    % table to contain the optimal matrices U
                               % rows indexed by N, columns indexed by m
TableMu(:,1) = ones(Nmax,1);   % \mu(1,N) = 1
for N = 1:Nmax
  TableMu(N,N) = 1;            % \mu(N,N) = 1
end
for N = 3:Nmax
  for m = 2:N-1
    [mu,B,U] = QMaxRelProjCst(m,N);
    TableMu(N,m) = mu;
    TableMuB{N,m} = B;
    TableMuU{N,m} = U;
  end
end

TableMuTilde = TableMu - 1./(1:N)'*(1:N);    % table of the reduced quantities \mu(m,N)-m/N

save('reproducible.mat','TableMu','TableMuB', 'TableMuU', 'TableMuTilde');


%% Verification that lambda(3,5) is strictly larger than mu(3,5)
% option 1: known value lambda(3,5) = (5+4*sqrt(2))/7

[(5+4*sqrt(2))/7  TableMu(5,3)]

%% option 2: calculate the projection constant of the space obtained in [1]
% note that one needs the software MinProj [2], which itself uses CVX [3]

a = sqrt((2*sqrt(2)-1)/7);
b = sqrt(1-a^2)/2;
u1 = sqrt((5-3*sqrt(2))/7);
c0 = sqrt((2+3*sqrt(2))/7)/2;

v1 = [a/u1; b/c0; b/c0; -b/c0; -b/c0];
v2 = [0; 0; 0; 1/sqrt(2)/c0; -1/sqrt(2)/c0];
v3 = [0; 1/sqrt(2)/c0; -1/sqrt(2)/c0; 0; 0];
V = [v1,v2,v3];

[MinProjCoor(V,inf) TableMu(5,3)]

%% option 3: lower-bound \lambda(3,5) using random choices of t and U
% warning: takes about 10min, output not guaranteed to exceed \mu(3,5) 

clear U;
tic;
lowb = 1;
for n = 1:25000000
   t = rand(5,1); t = t/norm(t);
   [U_aug,~,~] = svd(rand(5,5)); U = U_aug(:,1:3);
   lowb_aux = sum(sum(abs(diag(t)*U*U'*diag(t))));
   if lowb_aux > lowb
      lowb = lowb_aux; 
   end
end
[lowb TableMu(5,3)]
toc

%% option 4: lower-bound \lambda(3,5) by alternating maximizations
% see below where general lower bounds for \lambda(m,N) are obtained


%% Computation of the Real Maximal Relative Projection Constants \lambda(m,N)
% strickly speaking, we compute lower-bounds,
% but we believe them to be the true values.
% the function MaxRelProjCst(m,N) is defined in its separate file

Nmax = 10;
TableLambda = NaN(Nmax,Nmax);      % table to contain the values of \lambda
TableLambdaB = cell(Nmax,Nmax);    % table to contain the optimal Seidel matrices
TableLambdaT = cell(Nmax,Nmax);    % table to contain the optimal eigenvectors t
TableLambdaU = cell(Nmax,Nmax);    % table to contain the optimal matrices U
TableLambdaV = cell(Nmax,Nmax);    % table to contain the optimal spaces given by V
                                   % rows indexed by N, columns indexed by m
TableLambda(:,1) = ones(Nmax,1);   % \lambda(1,N) = 1
for N = 1:Nmax
  TableLambda(N,N) = 1;            % \lambda(N,N) = 1
end
for N = 3:Nmax
  for m = 2:N-1
    [lambda,B,t,U,V] = MaxRelProjCst(m,N);
    TableLambda(N,m) = lambda;
    TableLambdaB{N,m} = B;
    TableLambdaT{N,m} = t;
    TableLambdaU{N,m} = U;
    TableLambdaV{N,m} = V;
  end
end

TableLambdaTilde = TableLambda - 1./(1:N)'*(1:N);    % table of the reduced quantities \lambda(m,N)-m/N

save('reproducible.mat','TableLambda','TableLambdaB', 'TableLambdaT',...
    'TableLambdaU', 'TableLambdaV','TableLambdaTilde','-append');


%% Partial evidence that the lower bounds are the true values of \lambda(m,N)
% we apply the same alernating maximization strategy,
% but with different starting point than t^{(0)}=[1,...,1]/sqrt(N);
% this is embedded in the function MaxRelProjCst_randomized.
% We test the cases (m=4,N=7) and (m=5,N=8), 
% for which we should obtain 1.7398 and 1.8377, respectively.

nTests = 100;
MaxRelProjCst_randomized(4,7,nTests)
MaxRelProjCst_randomized(5,8,nTests)


%% Appendix A: verification of the formula obtained by trace duality
% for the projection constant of a given subspace
% note that one needs the software CVX [3] and the software MinProj [2], 
% which itself uses CVX

m =4; N = 7; 
V = randn(N,m);
pc1 = MinProjCoor(V,inf);  % the projection constant obtained with MinProj
cvx_begin
variable M(N,N);
variable C(m,m);
maximize trace(C)
subject to 
sum(max(abs(M'))) <= 1;
M'*V == V*C;
cvx_end
pc2 = cvx_optval;          % the projection constant obtained by trace duality
[pc1 pc2]                  % these two values should be similar

%% Appendix B: as described in the supplementary note, here is another way 
% to compute (lower bounds for) the (quasi)maximal relative projections constants 
% Consider first the real setting

Nmax = 10;                            % note: more values could be obtained 
nTest = 1000;                         % but one would have to increase nTest 
TableMu_Real = NaN(Nmax,Nmax);        % table to contain the values of \mu
TableLambda_Real = NaN(Nmax,Nmax);    % table to contain the values of \lambda
TableMu_Real(:,1) = ones(Nmax,1);     % \mu(1,N) = 1
TableLambda_Real(:,1) = ones(Nmax,1); % \lambda(1,N) = 1
for N = 1:Nmax
  TableMu_Real(N,N) = 1;              % \mu(N,N) = 1
  TableLambda_Real(N,N) = 1;          % \lambda(N,N) = 1
end
for N = 3:Nmax
  for m = 2:N-1
    mu = QMaxRelProjCst_Real(m,N,nTest);
    TableMu_Real(N,m) = mu;
    lambda = MaxRelProjCst_Real(m,N,nTest);
    TableLambda_Real(N,m) = lambda;
  end
end 

save('reproducible.mat','TableMu_Real', 'TableLambda_Real', '-append');

% the values just computed agree with the ones computed earlier,
% as certified by the fact that the following two quantities are small 
[max(max(abs(TableMu-TableMu_Real))) max(max(abs(TableLambda-TableLambda_Real)))]

%% Consider next the complex setting

Nmax = 10;                             
nTest = 1000;                          
TableMu_Complex = NaN(Nmax,Nmax);        % table to contain the values of \mu
TableLambda_Complex = NaN(Nmax,Nmax);    % table to contain the values of \lambda
TableMu_Complex(:,1) = ones(Nmax,1);     % \mu(1,N) = 1
TableLambda_Complex(:,1) = ones(Nmax,1); % \lambda(1,N) = 1
for N = 1:Nmax
  TableMu_Complex(N,N) = 1;              % \mu(N,N) = 1
  TableLambda_Complex(N,N) = 1;          % \lambda(N,N) = 1
end
for N = 3:Nmax
  for m = 2:N-1
    mu = QMaxRelProjCst_Complex(m,N,nTest);
    TableMu_Complex(N,m) = mu;
    lambda = MaxRelProjCst_Complex(m,N,nTest);
    TableLambda_Complex(N,m) = lambda;
  end
end
TableMuTilde_Complex = TableMu_Complex - 1./(1:N)'*(1:N);         % table of the reduced quantities \mu(m,N)-m/N
TableLambdaTilde_Complex = TableLambda_Complex - 1./(1:N)'*(1:N); % table of the reduced quantities \lambda(m,N)-m/N

save('reproducible.mat','TableMu_Complex', 'TableLambda_Complex',...
    'TableMuTilde_Complex', 'TableLambdaTilde_Complex', '-append');


%% References
% [1] B. L. Chalmers and G. Lewicki.
% Three-dimensional subspace of \ell_\infty^{(5)} with maximal projection
% constant.
% Journal of Functional Analysis 257.2 (2009): 553-592.
% 
% [2] S. Foucart
% Computation of minimal projections and extensions.
% Numerical Functional Analysis and Optimization 37.2 (2016): 159--185.
% 
% [3] CVX Research, Inc. 
% CVX: matlab software for disciplined convex programming, version 2.1. 
% http://cvxr.com/cvx, 2014.