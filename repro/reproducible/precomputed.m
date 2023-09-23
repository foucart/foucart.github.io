%% 
% Creates all real Seidel matrices up to switching and permuation equivalences.
% For a detailed explanation, see Section 3 of
% "On Maximal Relative Projection Constants"
% by S. Foucart and L. Skrzypek.

% We have used the software called nauty, available on 
% Brendan McKay's or Adolfo Piperno's webpages at
% http://cs.anu.edu.au/~bdm/nauty/
% http://pallini.di.uniroma1.it/
% to construct the adjacency matrices of all simple graphs on n = N-1 vertices,
% with N = 3, 4, ..., Nmax.
% Here is what has been done in nauty, with N = 5 for illustration: 
%
% Up to isomorphism, all simple graphs on n = 4 vertices are constructed
% by typing the following command in the terminal window: 
%
% geng 4 NonIsoGraphs04.g6
%
% so that the result is saved in a file with format g6.
% Then the graphs are converted into their (0,1)-adjancency matrices
% and the result is saved into a file with format txt by typing
%
% showg -Aq NonIsoGraphs04.g6 NonIsoGraphs04.txt
%
% Next, with a text editor, we manipulate the text file to replace all
% the occurence of n (i.e., 4 here) by blanks,
% so that the text file contains many nXn matrices separated by blank lines.
%
% Finally, we proceed in MATLAB to read the text file into a big matrix, 
% separate it into a bunch of (0,1)-adjacency matrices, 
% transfom them into nxn (-1,0,1)-adjacency matrices,
% and finish by creating NXN (-1,0,1)-adjacency matrices 
% by appending a first row/column of 1's (except the first entry) 
% to obtain all Seidel matrices up to equivalence.

% Written by Simon Foucart in June 2015, modified in February 2016
% Send comments to simon.foucart@centraliens.net

%%
Nmax = 10;
Seidel = cell(1,Nmax);

%%
N = 3;    % modify here when the size of Seidel matrices is changed
n = 2;    % modify here when the size of Seidel matrices is changed
I = eye(n);
J = ones(n,n);
BigMat = dlmread('NonIsoGraphs02.txt');    % modify here when the size of Seidel matrices is changed
NbMat = size(BigMat,1)/n;
for k = 1:NbMat
  temp = -2*BigMat((k-1)*n+1:k*n,:) + J - I;
  Seidel{N}(:,:,k) = [0, ones(1,n); ones(n,1), temp];
end

%%
N = 4;    % modify here when the size of Seidel matrices is changed
n = 3;    % modify here when the size of Seidel matrices is changed
I = eye(n);
J = ones(n,n);
BigMat = dlmread('NonIsoGraphs03.txt');    % modify here when the size of Seidel matrices is changed
NbMat = size(BigMat,1)/n;
for k = 1:NbMat
  temp = -2*BigMat((k-1)*n+1:k*n,:) + J - I;
  Seidel{N}(:,:,k) = [0, ones(1,n); ones(n,1), temp];
end

%%
N = 5;    % modify here when the size of Seidel matrices is changed
n = 4;    % modify here when the size of Seidel matrices is changed
I = eye(n);
J = ones(n,n);
BigMat = dlmread('NonIsoGraphs04.txt');    % modify here when the size of Seidel matrices is changed
NbMat = size(BigMat,1)/n;
for k = 1:NbMat
  temp = -2*BigMat((k-1)*n+1:k*n,:) + J - I;
  Seidel{N}(:,:,k) = [0, ones(1,n); ones(n,1), temp];
end

%%
N = 6;    % modify here when the size of Seidel matrices is changed
n = 5;    % modify here when the size of Seidel matrices is changed
I = eye(n);
J = ones(n,n);
BigMat = dlmread('NonIsoGraphs05.txt');    % modify here when the size of Seidel matrices is changed
NbMat = size(BigMat,1)/n;
for k = 1:NbMat
  temp = -2*BigMat((k-1)*n+1:k*n,:) + J - I;
  Seidel{N}(:,:,k) = [0, ones(1,n); ones(n,1), temp];
end

%%
N = 7;    % modify here when the size of Seidel matrices is changed
n = 6;    % modify here when the size of Seidel matrices is changed
I = eye(n);
J = ones(n,n);
BigMat = dlmread('NonIsoGraphs06.txt');    % modify here when the size of Seidel matrices is changed
NbMat = size(BigMat,1)/n;
for k = 1:NbMat
  temp = -2*BigMat((k-1)*n+1:k*n,:) + J - I;
  Seidel{N}(:,:,k) = [0, ones(1,n); ones(n,1), temp];
end

%%
N = 8;    % modify here when the size of Seidel matrices is changed
n = 7;    % modify here when the size of Seidel matrices is changed
I = eye(n);
J = ones(n,n);
BigMat = dlmread('NonIsoGraphs07.txt');    % modify here when the size of Seidel matrices is changed
NbMat = size(BigMat,1)/n;
for k = 1:NbMat
  temp = -2*BigMat((k-1)*n+1:k*n,:) + J - I;
  Seidel{N}(:,:,k) = [0, ones(1,n); ones(n,1), temp];
end

%%
N = 9;    % modify here when the size of Seidel matrices is changed
n = 8;    % modify here when the size of Seidel matrices is changed
I = eye(n);
J = ones(n,n);
BigMat = dlmread('NonIsoGraphs08.txt');    % modify here when the size of Seidel matrices is changed
NbMat = size(BigMat,1)/n;
for k = 1:NbMat
  temp = -2*BigMat((k-1)*n+1:k*n,:) + J - I;
  Seidel{N}(:,:,k) = [0, ones(1,n); ones(n,1), temp];
end

%% This is the Nmax used in this file
N = 10;    % modify here when the size of Seidel matrices is changed
n = 9;     % modify here when the size of Seidel matrices is changed
I = eye(n);
J = ones(n,n);
BigMat = dlmread('NonIsoGraphs09.txt');    % modify here when the size of Seidel matrices is changed
NbMat = size(BigMat,1)/n;
for k = 1:NbMat
  temp = -2*BigMat((k-1)*n+1:k*n,:) + J - I;
  Seidel{N}(:,:,k) = [0, ones(1,n); ones(n,1), temp];
end

%% save these precomputed Seidel matrices
save('precomputed.mat','Seidel');