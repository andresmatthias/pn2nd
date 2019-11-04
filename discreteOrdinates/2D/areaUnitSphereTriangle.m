function area = areaUnitSphereTriangle(A, B, C)
% AREAUNITSPHERETRIANGLE Compute the area of a geodesic triangle on the
%   unit sphere.
%   
%   The points A, B, C must be in same hemisphere; we get an error for 
%    A = [1, 0, 0]', B = [-1, 0, 0]', C = [0, 0, 1]'
%
%   Implementation based on l'Huillier theorem described in:
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

if checkUnitSphere(A) && checkUnitSphere(B) && checkUnitSphere(C)
   c = geodesicEdgeLengthsUnitSphere(A, B);
   a = geodesicEdgeLengthsUnitSphere(B, C);
   b = geodesicEdgeLengthsUnitSphere(A, C);  
   s = 1 / 2 * (a + b + c); % semi-perimeter 
   area = 4 * atan(sqrt(tan(s / 2) .* tan((s - a) / 2) .* tan((s - b) / 2) .* tan((s - c) / 2)));
end
end

function out = checkUnitSphere(A)
    out = all(abs(sqrt(sum(A.^2, 1)) -1) < 1e-15);
end