function [X, w] = gh_quad(N)
% GH_QUAD - N-point Gauss-Hermite quadrature
%
% SYNTAX: [X, w] = gh_quad(N)
%
% Returns the nodes X and weights w of the N-point
% Gauss-Hermite quadrature rule that integrates
% (w.r.t. standard Gaussian measure) exactly all
% polynomials up to degree 2N-1.
%
% INPUT
%   - N   the number of nodes
%
% OUTPUT
%   - X   the Gauss-Hermite nodes   (column vector)
%   - w   the Gauss-Hermite weights (column vector)

% Toni Karvonen, 2018
    
  % Use the Jacobi tridiagonal matrix to compute the
  % nodes and weights
  J = zeros(N, N);
  b = sqrt(1:N-1);
  J(N+1:N+1:end-1) = b;
  J(2:N+1:end-N) = b;
  [w,X] = eig(J);
  X = diag(X);
  w = (w(1,:).^2)';
  
  % Some trickery to output a completely symmetric node set
  Y = (abs(X(1:floor(N/2))) + flipud(X(ceil(N/2)+1:end)))/2;
  if mod(N,2) == 0;
    X = [-Y; flipud(Y)];
  else
    X = [-Y; 0; flipud(Y)];
  end

end
