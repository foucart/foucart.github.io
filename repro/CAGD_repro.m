%% Reproducible file accompanying the paper
% "Three topics in multivariate spline theory" by S. Foucart,
% specifically, its third section


%% Verification of the fact that 
% the symmetrized collocation matrix at the interior domain points ---
% i.e., the matrix C --- is positive definite up to d=17, but not for d=18

dmin = 3; dmax = 18;
lambdaC = zeros(1,dmax);
for d = dmin:dmax
  % definition of the matrix C
  C = zeros(nchoosek(d-1,2),nchoosek(d-1,2));
  for  i = 1:d-2
    for j = 1:d-1-i
      for mu = 1:d-2
        for nu = 1:d-1-mu
          C(nchoosek(d-1,2)-nchoosek(d-i,2)+j,nchoosek(d-1,2)-nchoosek(d-mu,2)+nu) ...
            = exp(d)*( nchoosek(d,i)*nchoosek(d-i,j)/d^d*mu^i*nu^j*(d-mu-nu)^(d-i-j) + ...
            nchoosek(d,mu)*nchoosek(d-mu,nu)/d^d*i^mu*j^nu*(d-i-j)^(d-mu-nu) );
        end
      end
    end
  end
  % determination of the smallest eignevalue of C
  lambdaC(d) = min(eig(C));
end
% display the list of smallest eigenvalues from dmin to dmax
[dmin:dmax;lambdaC(dmin:dmax)]

%% Validation of the reliability of the numerics
% by finding the minimal eigenvalue of the unsymmetrized collocation matrix 
% and comparing with the available theoretical results 

dmin = 3; dmax = 18;
lambdaB = zeros(1,dmax);
lambdaB_theor = zeros(1,dmax);
for d = dmin:dmax
  % definition of the matrix B
  B = zeros(nchoosek(d-1,2),nchoosek(d-1,2));
  for  i = 1:d-2
    for j = 1:d-1-i
      for mu = 1:d-2
        for nu = 1:d-1-mu
          B(nchoosek(d-1,2)-nchoosek(d-i,2)+j,nchoosek(d-1,2)-nchoosek(d-mu,2)+nu) ...
            = exp(d)*nchoosek(d,i)*nchoosek(d-i,j)*(1/d)^d*mu^i*nu^j*(d-mu-nu)^(d-i-j);
        end
      end
    end
  end
  % determination of the smallest eignevalue of B, experimentally and
  % theoretically
  lambdaB(d) = min(real(eig(B)));
  lambdaB_theor(d) = exp(d)*factorial(d)/d^d;
end
% display the lists of experimental and theoretical smallest eigenvalues
[dmin:dmax; lambdaB(dmin:dmax); lambdaB_theor(dmin:dmax)]


%% Verification of the fact that, for the trivariate case, 
% the symmetrized collocation matrix at the interior domain points ---
% i.e., the matrix C --- is positive definite up to d=16, but not for d=17

dmin = 4; dmax = 17;
lambdaC3 = zeros(1,dmax);
for d = dmin:dmax
  % definition of the matrix C
  C3 = zeros(nchoosek(d-1,3),nchoosek(d-1,3));
  for  i = 1:d-3
    for j = 1:d-2-i
      for k = 1:d-1-i-j
        for mu = 1:d-3
          for nu = 1:d-2-mu
            for kappa = 1:d-1-mu-nu
              C3(nchoosek(d-1,3)-nchoosek(d-i,3)+nchoosek(d-i-1,2)-nchoosek(d-i-j,2)+k,...
                nchoosek(d-1,3)-nchoosek(d-mu,3)+nchoosek(d-mu-1,2)-nchoosek(d-mu-nu,2)+kappa) ...
                = exp(d)*( nchoosek(d,i)*nchoosek(d-i,j)*nchoosek(d-i-j,k)/d^d...
                *mu^i*nu^j*kappa^k*(d-mu-nu-kappa)^(d-i-j-k) + ...
                nchoosek(d,mu)*nchoosek(d-mu,nu)*nchoosek(d-mu-nu,kappa)/d^d...
                *i^mu*j^nu*k^(kappa)*(d-i-j-k)^(d-mu-nu-kappa) );
            end
          end
        end
      end
    end
  end
  % determination of the smallest eignevalue of C
  lambdaC3(d) = min(eig(C3));
end
% display the list of smallest eigenvalues from dmin to dmax
[dmin:dmax;lambdaC3(dmin:dmax)]

%% Validation of the reliability of the numerics
% by finding the minimal eigenvalue of the unsymmetrized collocation matrix 
% and comparing with the available theoretical results 

dmin = 4; dmax = 17;
lambdaB3 = zeros(1,dmax);
for d = dmin:dmax
  % definition of the matrix B
  B3 = zeros(nchoosek(d-1,3),nchoosek(d-1,3));
  for  i = 1:d-3
    for j = 1:d-2-i
      for k = 1:d-1-i-j
        for mu = 1:d-3
          for nu = 1:d-2-mu
            for kappa = 1:d-1-mu-nu
              B3(nchoosek(d-1,3)-nchoosek(d-i,3)+nchoosek(d-i-1,2)-nchoosek(d-i-j,2)+k,...
                nchoosek(d-1,3)-nchoosek(d-mu,3)+nchoosek(d-mu-1,2)-nchoosek(d-mu-nu,2)+kappa) ...
                = exp(d)*nchoosek(d,i)*nchoosek(d-i,j)*nchoosek(d-i-j,k)/d^d...
                *mu^i*nu^j*kappa^k*(d-mu-nu-kappa)^(d-i-j-k);
            end
          end
        end
      end
    end
  end
  % determination of the smallest eignevalue of B, experimentally and
  % theoretically
  lambdaB3(d) = min(real(eig(B3)));
  lambdaB3_theor(d) = exp(d)*factorial(d)/d^d;
end
% display the lists of experimental and theoretical smallest eigenvalues
[dmin:dmax; lambdaB3(dmin:dmax); lambdaB3_theor(dmin:dmax)]


%% Verification of the fact that, in the quadrivariate case,
% the symmetrized collocation matrix at the interior domain points ---
% i.e., the matrix C --- is positive definite up to d=14, but not for d=15

dmin = 5; dmax = 15;
lambdaC4 = zeros(1,dmax);
for d = dmin:dmax
  % definition of the matrix C
  C4 = zeros(nchoosek(d-1,4),nchoosek(d-1,4));
  for  i = 1:d-4
    for j = 1:d-3-i
      for k = 1:d-2-i-j
        for l = 1:d-1-i-j-k
          for mu = 1:d-4
            for nu = 1:d-3-mu
              for kappa = 1:d-2-mu-nu
                for lambda = 1:d-1-mu-nu-kappa
                  idxrow=nchoosek(d-1,4)-nchoosek(d-i,4)+nchoosek(d-i-1,3)-nchoosek(d-i-j,3)...
                    +nchoosek(d-i-j-1,2)-nchoosek(d-i-j-k,2)+l;
                  idxcol=nchoosek(d-1,4)-nchoosek(d-mu,4)+nchoosek(d-mu-1,3)-nchoosek(d-mu-nu,3)...
                    +nchoosek(d-mu-nu-1,2)-nchoosek(d-mu-nu-kappa,2)+lambda;
                  C4(idxrow,idxcol) ...
                    = exp(d)*( nchoosek(d,i)*nchoosek(d-i,j)*nchoosek(d-i-j,k)*nchoosek(d-i-j-k,l)/d^d...
                    *mu^i*nu^j*kappa^k*lambda^l*(d-mu-nu-kappa-lambda)^(d-i-j-k-l) + ...
                    nchoosek(d,mu)*nchoosek(d-mu,nu)*nchoosek(d-mu-nu,kappa)*nchoosek(d-mu-nu-kappa,lambda)/d^d...
                    *i^mu*j^nu*k^(kappa)*l^lambda*(d-i-j-k-l)^(d-mu-nu-kappa-lambda) );
                end
              end
            end
          end
        end
      end
    end
  end
  % determination of the smallest eignevalue of C
  lambdaC4(d) = min(eig(C4));
end
% display the list of smallest eigenvalues from dmin to dmax
[dmin:dmax;lambdaC4(dmin:dmax)]

%% Validation of the reliability of the numerics
% by finding the minimal eigenvalue of the unsymmetrized collocation matrix 
% and comparing with the available theoretical results 

dmin = 5; dmax = 15;
lambdaB4 = zeros(1,dmax);
for d = dmin:dmax
  % definition of the matrix B
  B4 = zeros(nchoosek(d-1,4),nchoosek(d-1,4));
  for  i = 1:d-4
    for j = 1:d-3-i
      for k = 1:d-2-i-j
        for l = 1:d-1-i-j-k
          for mu = 1:d-4
            for nu = 1:d-3-mu
              for kappa = 1:d-2-mu-nu
                for lambda = 1:d-1-mu-nu-kappa
                  idxrow=nchoosek(d-1,4)-nchoosek(d-i,4)+nchoosek(d-i-1,3)-nchoosek(d-i-j,3)...
                    +nchoosek(d-i-j-1,2)-nchoosek(d-i-j-k,2)+l;
                  idxcol=nchoosek(d-1,4)-nchoosek(d-mu,4)+nchoosek(d-mu-1,3)-nchoosek(d-mu-nu,3)...
                    +nchoosek(d-mu-nu-1,2)-nchoosek(d-mu-nu-kappa,2)+lambda;
                  B4(idxrow,idxcol) ...
                    = exp(d)*nchoosek(d,i)*nchoosek(d-i,j)*nchoosek(d-i-j,k)*nchoosek(d-i-j-k,l)/d^d...
                    *mu^i*nu^j*kappa^k*lambda^l*(d-mu-nu-kappa-lambda)^(d-i-j-k-l) ;
                end
              end
            end
          end
        end
      end
    end
  end
  % determination of the smallest eignevalue of B, experimentally and
  % theoretically
  lambdaB4(d) = min(real(eig(B4)));
  lambdaB4_theor(d) = exp(d)*factorial(d)/d^d;
end
% display the lists of experimental and theoretical smallest eigenvalues
[dmin:dmax; lambdaB4(dmin:dmax); lambdaB4_theor(dmin:dmax)]