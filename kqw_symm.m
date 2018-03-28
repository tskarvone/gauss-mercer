function w = kqw_symm(X, k, kmean)
% KQW_SYMM - weights of kernel quadrature
%
% SYNTAX: w = kqw_symm(X, k, kmean)
%
% Computes the weights w of a (univariate) kernel quadrature rule 
% for the kernel k with the kernel mean function kmean at a 
% symmetric 1D nodeset X (i.e. both x and -x are in X).
%
% The weight are computed using the algorithm presented in 
%
%   T. Karvonen and S. Särkkä (2018). Fully symmetric kernel quadrature. 
%   SIAM Journal on Scientific Computing, 40(2):A697–A720.
%   https://doi.org/10.1137/17M1121779
%
% This makes the computation (slightly) more numerically stable and
% ensures that the weights for -x and x are equal.
%
% INPUT
%   - X           a vector of symmetric 1D nodes
%   - k           kernel, given as a bivariate function k(x,y)
%   - kmean       kernel mean, given as kmean(x)
%
% OUTPUT
%   - w           column vector of the weights

% Toni Karvonen, 2018
  
  N = length(X);
  n = floor(N/2);
  XX = X(1:n);
  
  if mod(N,2) == 0
    KK = zeros(n,n);
    kmv = zeros(n,1);
    for i = 1:n
      kmv(i) = kmean(XX(i));
      for j = 1:n
        KK(i,j) = k(XX(i),XX(j)) + k(-XX(i),XX(j));
      end
    end
    w = KK\kmv;
    w = [w; flipud(w)];
  else
    KK = zeros(n+1,n+1);
    kmv = zeros(n+1,1);
    for i = 1:n
      kmv(i) = kmean(XX(i));
      for j = 1:n
        KK(i,j) = k(XX(i),XX(j)) + k(-XX(i),XX(j));
      end
    end
    
    kmv(end) = kmean(0);
    for i = 1:n
      KK(n+1,i) = 2 * k(0,XX(i));
      KK(i,n+1) = k(XX(i),0);
    end
    KK(n+1,n+1) = k(0,0);
    w = KK\kmv;
    w = [w; flipud(w(1:end-1))];
  end

end
