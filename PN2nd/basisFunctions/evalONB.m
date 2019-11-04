function [basisAtPoints] = evalONB(pointsOnSphereMuPhi, N, spatialDimension)
% EVALONB Evaluate orthonormal basis functions (real spherical harmonics)
%   up to degree N at given points on sphere.
%
%   Depending on spatialDimension neglect certain basis functions due to
%   symmetry assumptions.
%
%   Parametrization of points on sphere:
%       Omega = [cos(phi) * sqrt(1 - mu^2);
%                sin(phi) * sqrt(1 - mu^2);
%                mu]
% 
%   For a demonstration on how to use this function,
%   see also BASISFUNCTIONSTEST, FLUXMATRIXTEST
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

basisAtPoints = evalRealSphericalHarmonics(pointsOnSphereMuPhi, N);
fullSpatialDimension = 3;
nMomentsFull = getNumberOfBasisFunctions(N, fullSpatialDimension);
linIdx = 0 : nMomentsFull - 1;
[l, m] = linearIdx2DegOrder(linIdx, fullSpatialDimension);
idxReduced = getReducedIdx(l, m, spatialDimension); 
basisAtPoints = basisAtPoints(idxReduced, :); 
end