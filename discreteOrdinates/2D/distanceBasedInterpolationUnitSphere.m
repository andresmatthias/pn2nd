function [idxOrdinates, weights] = distanceBasedInterpolationUnitSphere(pointsOnSphereCartesian, queryPointsCartesian)
% DISTANCEBASEDINTERPOLATIONUNITSPHERE Get weights for interpolating a
%   function on the unit sphere. The weights are computed via the geodesic
%   distances of the query point to the reference points.
%   
%   Implementation based on:
%       Interpolation on spherical geodesic grids: A comparative study; 
%       Maria Francesca Carfora; 2007; Journal of Computational and Applied 
%       Mathematics; Volume 210;
%       url: https://www.sciencedirect.com/science/article/pii/S0377042706006522
% 
%   For a demonstration on how to use this function,
%   see also DISCRETEORDINATES2DTEST
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

idxOrdinates = zeros(3, size(queryPointsCartesian, 2));
weights = zeros(3, size(queryPointsCartesian, 2));
for k = 1 : size(weights, 2)
    arcLength = geodesicEdgeLengthsUnitSphere(pointsOnSphereCartesian, queryPointsCartesian(:, k));
    [m, idx] = mink(arcLength, 3);
    if m(1) < 1e-14
       w = zeros(3, 1);
       w(1) = 1;
    else
        w = 1 ./ arcLength(idx); 
    end
    weights(:, k) = w / sum(w);
    idxOrdinates(:, k) = idx;
end
end

