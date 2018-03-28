function Y = hermite_probn(n, X)
% HERMITE_PROBN - normalised probabilists' Hermite polynomials
%
% SYNTAX: Y = hermite_probn(n,X)
%
% Returns element-wise evaluations Y of the nth normalised 
% probabilist's Hermite polynomial at points X.
%
% INPUT
%   - n     degree of the Hermite polynomial
%   - X     the evaluation points
%
% OUTPUT
%   - Y     the nth normalised probabilist's Hermite
%           polynomial evaluated at X

% Toni  Karvonen, 2018
  
  N = length(X);
  
  if n == 0
    Y = ones(size(X));
    return
  end
  
  if n == 1
    Y = X;
    return
  end
  
  % Use the recurrence relation
  Yp = ones(size(X));
  Y = X;
  for m = 2:n
    tmp = Y;
    Y = X.*Y - (m-1)*Yp;
    Yp = tmp;
  end
  Y = Y/sqrt(factorial(n));
  
end
