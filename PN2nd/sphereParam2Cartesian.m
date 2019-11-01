function Omega = sphereParam2Cartesian(mu, phi)
% SPHEREPARAM2CARTESIAN Represent a point on the unit sphere, given in 
%   angular parametrization (height-azimuthal), in Cartesian coordinates.
%
%   For a demonstration on how to use this function,
%   see also MISCTEST, FLUXMATRIXTEST
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

mu = mu(:)';
phi = phi(:)';
Omega = [sqrt(1 - mu.^2) .* cos(phi);
         sqrt(1 - mu.^2) .* sin(phi);
         mu];
end