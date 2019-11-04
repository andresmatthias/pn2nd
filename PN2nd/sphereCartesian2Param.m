function [mu, phi] = sphereCartesian2Param(pointsOnSphere)
% SPHERECARTESIAN2PARAM Represent a point on the unit sphere, given in 
%   Cartesian coordinates, in angular parametrization (height-azimuthal).
%
%   For a demonstration on how to use this function,
%   see also MISCTEST
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

mu = pointsOnSphere(3, :);
phi = atan2(pointsOnSphere(2, :), pointsOnSphere(1, :));
end