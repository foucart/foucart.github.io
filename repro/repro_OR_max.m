%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reproducible MATLAB file accompanying the paper
% LEARNING THE MAXIMUM OF A HOLDER FUNCTION FROM INEXACT DATA
% by S. Foucart
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Written by Simon Foucart in February 2025
% Send commments to simon.foucart@centraliens.net

clear all; clc;

%% Create a Lipschitz function and produce inexact point values

% the Lipschitz function f
a = 0.8; b = 0.2;
f = @(x) abs( a/5*sin(5*x) + b*(x.^4-1.5*x.^2) )/2 ...
       + abs( a/2*cos(2*x) + b/50*sin(50*x) )/2;

% a fine grid for the interval [-1,1]
grid_size = 1000; 
grid = linspace(-1,1,grid_size)';

% the datasites x_m --- for later convenience, they are taken on the grid
M = 9; % number of datasites
aux = sort(randi(grid_size,M,1));
X = grid(aux);  

% the inexact point values y_m = f(x_m) + e_m
Eps = 0.04*ones(M,1); % here, the bounds epsilon_m on |e_m| are all equal
E = Eps.*(2*rand(M,1)-1);
Y = f(X) + E;

% visualization of f and of the data
figure(1)
plot(grid,f(grid),'b',X,Y,'ko')
ylim([-0.3 0.5])
hold on
figure(2)
plot(grid,f(grid),'b',X,Y,'ko')
ylim([-0.3 0.5])
hold on

%% Produce and visualize the lower function l, the upper function u, and
% the locally optimal estimation (l+u)/2 for full approximation
% (the definitions of l and u are at the end of the file)

l = lb(X,Y,Eps,grid);
u = ub(X,Y,Eps,grid);
loc = (l+u)/2;
figure(1)
plot(grid,l,'k:',grid,u,'k:',grid,loc,'r--')

%% Produce and visualize the locally optimal estimation for the maximum,
% as well as the globally optimal estimation for the maximum

loc_max = ( max(l) + max(u) )/2;
U = ub(X,zeros(M,1),zeros(M,1),grid);
glo_max = max(Y) + max(U)/2;
figure(2)
plot(grid,loc_max*ones(grid_size,1),'r--',grid,glo_max*ones(grid_size,1),'g-.')

%% Definitions of the lower and upper functions

function res = lb(X,Y,Eps,grid)

M = length(X);
N = length(grid);
aux = zeros(M,N);
for m = 1:M
    for n=1:N
        aux(m,n) = Y(m) - Eps(m) - abs(grid(n)-X(m));
    end
end
res = zeros(N,1);
for n = 1:N
    res(n) = max(aux(:,n));
end

end 

function res = ub(X,Y,Eps,grid)

M = length(X);
N = length(grid);
aux = zeros(M,N);
for m = 1:M
    for n=1:N
        aux(m,n) = Y(m) + Eps(m) + abs(grid(n)-X(m));
    end
end
res = zeros(N,1);
for n = 1:N
    res(n) = min(aux(:,n));
end

end