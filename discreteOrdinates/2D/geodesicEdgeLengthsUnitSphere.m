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
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
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