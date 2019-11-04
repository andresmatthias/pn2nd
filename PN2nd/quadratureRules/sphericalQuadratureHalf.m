function [weights, mu, phi] = sphericalQuadratureHalf(outerNormalVector, maxExactDegree)
% SPHERICALQUADRATUREHALF Get a quadrature rule for functions defined on
%   the half of the unit sphere w.r.t. an outer normal vector n, i.e.,
%   (n * Omega_z < 0).
%
%   The quadrature nodes are given in the following parametrization of
%   the unit sphere:
%   Omega = [x; y; z] = [cos(phi) * sqrt(1 - mu^2);
%                        sin(phi) * sqrt(1 - mu^2);
%                        mu]
%
%   It is exact for polynomials on the sphere with degree <= maxExactDegree
%   (w.r.t. to Cartesian coordinates, like the real spherical harmonics, 
%   e.g., S_1^{-1} ~ y, S_2^{-2} ~ xy, S_2^2    ~ x^2 - y^2) 
%
%   For a demonstration on how to use this function,
%   see also SPHERICALQUADRATURETEST
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

% Rotate the nodes of the quadrature rule for the upper half with the same
% transformation that rotates -e3 = [0; 0; -1] to the outer normal vector.
addpath('../')
[weights, mu0, phi0] = sphericalQuadratureUpperHalf(maxExactDegree);
Omega0 = sphereParam2Cartesian(mu0, phi0);
rotationMatrix = getGeomRotationa2b([0; 0; -1], outerNormalVector);
Omega = rotationMatrix * Omega0;
[mu, phi] = sphereCartesian2Param(Omega);
end