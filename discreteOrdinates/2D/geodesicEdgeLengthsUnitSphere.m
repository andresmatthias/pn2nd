function arcLength = geodesicEdgeLengthsUnitSphere(A, B)
% GEODESICEDGELENGTHUNITSPHERE Compute the length of the geodesic edge
%   between two points on the unit sphere.
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

    if onUnitSphereCheck(A) && onUnitSphereCheck(B)
        arcLength = 2 * asin(sqrt(sum((A - B).^2, 1)) / 2);
    else
       error('Points not normalized, thus not on unit sphere!') 
    end
end

function out = onUnitSphereCheck(A)
    out = all(abs(sum(A.^2, 1) - 1) < 1e-14);
end