%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reproducible file accompanying the paper
%   FACILITATING OWL NORM MINIMIZATIONS
% by S. Foucart
% Written by S. Foucart in December 2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;
% CVX [1] is required to execute this reproducible
cvx_quiet true
cvx_solver gurobi     % remove if you do not have the gurobi license 

%% Verification of the linear characterizations of the OWL norm
% i.e., numerical validation of Theorem 1

n = 100;
x = randn(n,1);
w = sort(rand(n,1),'descend');

% OWL norm of x via its defining expression
x_star = sort(abs(x),'descend');
OWL = w'*x_star;

% OWL norm of x via the linear maximization
tic;
cvx_begin
variable S(n,n)
maximize sum(w.*(S*abs(x)))
subject to
S >= 0;
sum(S,1) == ones(1,n);
sum(S,2) == ones(n,1);
cvx_end
OWL_max = cvx_optval;
t_max = toc;

% OWL norm of x via the linear minimization
tic;
cvx_begin
variable a(n)
variable b(n)
minimize sum(a)+sum(b)
subject to
repmat(a,1,n) + repmat(b',n,1) >= +w*x';
repmat(a,1,n) + repmat(b',n,1) >= -w*x';
cvx_end
OWL_min = cvx_optval;
t_min = toc;

% comparison of the three outpouts
sprintf(strcat(...
    'The defining expression yields a value of', 32, num2str(OWL), 10,...
    'The linear maximization yields a value of', 32, num2str(OWL_max),...
    ', computed in', 32, num2str(t_max), ' seconds', 10,...
    'The linear minimization yields a value of', 32, num2str(OWL_min),...
    ', computed in', 32, num2str(t_min), ' seconds'))


%% Verification of the linear characterizations of the Ky-Fan norm
% i.e., numerical validation of Theorem 2
% note: these linear programs are indeed faster than their OWL counterparts

k = floor(n/2);
w = [ones(k,1); zeros(n-k,1)];

% Ky-Fan norm of x via its defining expression
KF = sum(x_star(1:k));

% Ky-Fan norm of x via the linear maximization
tic;
cvx_begin
variable u(n)
variable v(n)
maximize sum(v.*x)
subject to
sum(u) <= k;
u - v >= 0;
u + v >= 0;
u <= 1;
cvx_end
KF_max = cvx_optval;
t_max = toc;

% Ky-Fan norm of x via the linear minimization
tic;
cvx_begin
variable a(n)
variable b(n)
variable alp(n)
variable bet
minimize sum(alp)+k*bet
subject to
a + b == x;
alp - a >= 0;
alp + a >= 0;
bet - b >= 0;
bet + b >= 0;
cvx_end
KF_min = cvx_optval;
t_min = toc;

% comparison of the three outpouts
sprintf(strcat(...
    'The defining expression yields a value of', 32, num2str(KF), 10,...
    'The linear maximization yields a value of', 32, num2str(KF_max),...
    ', computed in', 32, num2str(t_max), ' seconds', 10,...
    'The linear minimization yields a value of', 32, num2str(KF_min),...
    ', computed in', 32, num2str(t_min), ' seconds'))


%% Illustration: solving an OWL norm minization under equality constraints
% If the last two columns of the measurement matrix are identical,
% then the last two entries of the solution are forced to be equal 

n = 200;                             % ambient dimension
m = 100;                             % number of measurements
s = 10;                              % sparsity level
% construct a measurement matrix with last two columns being identical
A_aux = randn(m,n-1);
A = [A_aux A_aux(:,n-1)];
% construct a sparse vector with last two entries being equal
supp_aux = randperm(n-2);  supp = sort(supp_aux(1:s-1));
x = zeros(n,1);
x(supp) = randn(s-1,1); 
x(n) = 1/2; x(n-1) = 1/2;
% produce the compressive measurement vector
y = A*x;                             

% attempt to recover x from y by L1 norm minimization
cvx_begin
variable xL1(n);
variable c(n);
minimize sum(c)
subject to 
A*xL1 == y;
c + xL1 >= 0;
c - xL1 >= 0;
cvx_end

% attemp to recover x from y by OWL norm minimization
w = sort(rand(n,1),'descend');
cvx_begin
variable xOWL(n);
variable a(n);
variable b(n);
minimize sum(a)+sum(b)
subject to
A*xOWL == y;
repmat(a,1,n) + repmat(b',n,1) >= +w*xOWL';
repmat(a,1,n) + repmat(b',n,1) >= -w*xOWL';
cvx_end

sprintf(strcat('Recovery by L1 norm minimization is unsuccesful:', 32,...
    'the L1 error is', 32, num2str(norm(x-xL1,1)),10,...
    'Recovery by OWL norm minimization is succesful:', 32, 32,...
    'the L1 error is', 32, num2str(norm(x-xOWL,1))))

%% Substantiation of the linear characterizations of the dual OWL norm
% i.e., numerical validation of Theorem 3

n = 100;
x = randn(n,1);
w = sort(rand(n,1),'descend');

% dual OWL norm of x via the linear maximization
tic;
cvx_begin
variable z(n)
variable a(n)
variable b(n)
maximize sum(x.*z)
subject to 
sum(a)+sum(b) <= 1;
repmat(a,1,n) + repmat(b',n,1) >= +w*z';
repmat(a,1,n) + repmat(b',n,1) >= -w*z';
cvx_end
DOWL_max = cvx_optval;
t_max = toc;

% dual OWL norm of x via the linear maximization
tic;
cvx_begin
variable c
variable U(n,n)
variable V(n,n)
minimize c
subject to 
U >= 0;
V >= 0;
(U-V)*w == x;
sum(U+V,1) == c*ones(1,n);
sum(U+V,2) == c*ones(n,1);
cvx_end
DOWL_min = cvx_optval;
t_min = toc;

% comparison of the two outpouts
sprintf(strcat(...
    'The linear maximization yields a value of', 32, num2str(DOWL_max),...
    ', computed in', 32, num2str(t_max), ' seconds', 10,...
    'The linear minimization yields a value of', 32, num2str(DOWL_min),...
    ', computed in', 32, num2str(t_min), ' seconds'))

%% Verification of the linear characterizations of the dual Ky-Fan norm
% i.e., numerical validation of Theorem 4

k = floor(n/2);
w = [ones(k,1); zeros(n-k,1)];

% Dual OWL norm of x via its abridged expression                            
DKF = max(norm(x,Inf),norm(x,1)/k);

% Dual OWL norm of x via the linear maximization
tic;
cvx_begin
variable z(n)
variable a(n)
variable b(n)
variable alp(n)
variable bet
maximize sum(x.*z)
subject to 
sum(alp)+k*bet <= 1;
a + b == z;
alp - a >= 0;
alp + a >= 0;
bet - b >= 0;
bet + b >= 0;
cvx_end
DKF_max = cvx_optval;
t_max = toc;

% Dual OWL norm of x via the linear minimization
tic;
cvx_begin
variable c
variable alp 
variable bet(n)
minimize c
subject to 
alp <= c;
sum(bet) <= k*c;
alp - x >= 0;
alp + x >= 0;
bet - x >= 0;
bet + x >= 0;
cvx_end
DKF_min = cvx_optval;
t_min = toc;

% comparison of the three outpouts
sprintf(strcat(...
    'The abridged expression yields a value of', 32, num2str(DKF), 10,...
    'The linear maximization yields a value of', 32, num2str(DKF_max),...
    ', computed in', 32, num2str(t_max), ' seconds', 10,...
    'The linear minimization yields a value of', 32, num2str(DKF_min),...
    ', computed in', 32, num2str(t_min), ' seconds'))


%% References

% 1. CVX Research, Inc., 
% "CVX: MATLAB software for disciplined convex programming"
% version 2.1, 2014, http://cvxr.com/cvx.
