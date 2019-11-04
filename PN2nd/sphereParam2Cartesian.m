function Omega = sphereParam2Cartesian(mu, phi)
% SPHEREPARAM2CARTESIAN Represent a point on the unit sphere, given in 
%   angular parametrization (height-azimuthal), in Cartesian coordinates.
%
%   For a demonstration on how to use this function,
%   see also MISCTEST, FLUXMATRIXTEST
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

mu = mu(:)';
phi = phi(:)';
Omega = [sqrt(1 - mu.^2) .* cos(phi);
         sqrt(1 - mu.^2) .* sin(phi);
         mu];
end