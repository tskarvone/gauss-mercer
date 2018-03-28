function w = kqw(X, k, kmean)
% KQW - weights of kernel quadrature
%
% SYNTAX: w = kqw(X, k, kmean)
%
% Computes the weights w of a (univariate) kernel quadrature rule 
% for the kernel k with the kernel mean function kmean at the nodes X.
%
% INPUT
%   - X           a vector of 1D nodes
%   - k           kernel, given as a bivariate function k(x,y)
%   - kmean       kernel mean, given as kmean(x)
%
% OUTPUT
%   - w           column vector of the weights

% Toni Karvonen, 2017-2018
  
  N = length(X);
  K = zeros(N,N);
  kmv = zeros(N,1);
  for i = 1:N
    kmv(i) = kmean(X(i));
    for j = 1:N
      K(i,j) = k(X(i),X(j));
    end
  end
  w = K\kmv;

end
