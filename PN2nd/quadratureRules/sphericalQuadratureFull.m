function [weights, mu, phi] = sphericalQuadratureFull(maxExactDegree)
% SPHERICALQUADRATUREFULL Get a quadrature rule for surface integrals on 
%   the unit sphere. Polynomials up to maxExactDegree will be integrated
%   exactly up to numerical precision. 
%
%   The quadrature nodes are given in the following parametrization of
%   the unit sphere:
%   Omega = [x; y; z] = [cos(phi) * sqrt(1 - mu^2);
%                        sin(phi) * sqrt(1 - mu^2);
%                        mu]
%                     = [cos(phi) * sin(theta);
%                        sin(phi) * sin(theta);
%                        cos(theta)]
%
%   It is exact for polynomials with degree <= maxExactDegree, 
%   w.r.t. to Cartesian coordinates, like the real spherical harmonics, 
%   e.g., S_1^{-1} ~ y, S_2^{-2} ~ xy, S_2^2 ~ x^2 - y^2. 
% 
%   Based on the representation as trigonometric polynomial of degree n:
%   f = c0 + c1*z + ... + cn*z^n; restrict to unit circle with 
%   z = e^{i*k*z} for -n <= k<= n, and obtain trigonometric polynomial
%   f = a0 + sum_{k=1}^n a_k * cos(k*phi) + b_k * sin(k*phi). Then combine
%   quadrature rules for trigonometric polynomials of degree n, depending 
%   on phi and theta.
% 
%   Surface integral in theta, phi:
%   int_{0}^{pi} \int_{0}^{2*pi} f(theta, phi) * sin(theta) dphi dtheta
%   If, for fixed phi, f is a trigonometric polynomial in theta with 
%   degree <= n, then f(theta, phi) * sin(theta) is of degree <= n + 1.
%
%   For a demonstration on how to use this function,
%   see also SPHERICALQUADRATURETEST, FLUXMATRIXTEST, BASISFUNCTIONSTEST
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

% quadrature has order p => all polynomials with degree <= p-1 will be 
% integrated exactly
% 

tw = trigauss(maxExactDegree, 0, 2 * pi);
phi = tw(:, 1);
wphi = tw(:, 2);

% Attention: add 1 to maxExactDegree in theta due to functional
% determinant (see above)!
tw = trigauss(maxExactDegree + 1, 0, pi);
theta = tw(:, 1);
wtheta = tw(:, 2);
mu = cos(theta);
wmu = abs(wtheta .* sin(theta)); % * functional determinant

[phi, mu] = meshgrid(phi, mu);
[wphi, wmu] = meshgrid(wphi, wmu);
weights = wmu(:) .* wphi(:);
mu = mu(:)';
phi = phi(:)';
weights = weights(:)';

end

