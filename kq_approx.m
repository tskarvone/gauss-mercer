function [X, w] = kq_approx(l,a,N)
% KQ_APPROX - eigendecomposition to the kernel quadrature
%             weights when the kernel and measure are Gaussian
%
% SYNTAX: [X, w] = kq_approx(l,a,N)
%
% Returns N scaled Gauss-Hermite nodes X and approximate kernel
% quadrature weights w for the Gaussian kernel with length-scale
% l and for Mercer eigenfunctions orthogonal w.r.t. the Gaussian
% measure with variance 1/(2a^2). The integration measure is the
% standard Gaussian.
%
% See Theorem 2.3 for details.
%
% INPUT:
%   - l   length-scale of the Gaussian kernel
%   - a   global scale parameter; the Mercer eigenfunctions are
%         orthogonal w.r.t. the Gaussian measure with variance 1/(2a^2)
%         RECOMMENDED VALUE: a = 1/sqrt(2)
%   - N   the number of nodes
%
% OUTPUT
%   - X   the scaled Gauss-Hermite nodes                     (column vector)
%   - w   the approximate Gaussian kernel quadrature weights (column vector)

% Toni Karvonen, 2018
  
  % Compute the necessary constants
  ep = 1/(sqrt(2)*l);
  b = (1+(2*ep/a)^2)^(1/4);
  d2 = (a^2/2)*(b^2-1);
  
  % Construct the nodes by scaling the Gauss-Hermite nodes
  [Xg, wg] = gh_quad(N);
  X = Xg/(sqrt(2)*a*b);

  % Use the main theorem to compute the approximate weights
  ev = exp(d2*X.^2);
  Nf = floor((N-1)/2)+1;
  Vh = zeros(N,Nf);
  for m = 1:Nf
    Vh(:,m) = hermite_probn(2*(m-1),Xg);
  end
  
  Ns = (0:Nf-1)';
  pmuh = sqrt(factorial(2*Ns))./(2.^Ns.*factorial(Ns)) .* (2*a^2*b^2/(1+2*d2) - 1).^Ns;
  w =  sqrt(1/(1+2*d2)) * wg .* ev .* (Vh*pmuh);
  
end
