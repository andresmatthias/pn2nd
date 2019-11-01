function [mu, phi] = sphereCartesian2Param(pointsOnSphere)
% SPHERECARTESIAN2PARAM Represent a point on the unit sphere, given in 
%   Cartesian coordinates, in angular parametrization (height-azimuthal).
%
%   For a demonstration on how to use this function,
%   see also MISCTEST
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

mu = pointsOnSphere(3, :);
phi = atan2(pointsOnSphere(2, :), pointsOnSphere(1, :));
end