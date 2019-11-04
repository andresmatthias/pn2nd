function [] = generateTestCase6GeoFile()
% GENERATETESTCASE6GEOFILE Generate the .geo file (Gmsh) for test case 6.
% 
%   For the use in the second-order formulation of the PN equations, we
%   need to assign a new boundary ID for every occuring normal vector 
%   within the same part of the boundary.
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

addpath('../../tools/')
meshName = './meshTestCase6';
%%
radiusOuter = 1.5;
radiusInner = 1;
geoGenParams.nTop = 5;
geoGenParams.nLeft = 5;
geoGenParams.nInnerOuter = 20;
phi = linspace(3/2*pi, 2 * pi, geoGenParams.nInnerOuter);
m = [0, radiusOuter];


polyOuter = m + radiusOuter * [cos(phi)', sin(phi)'];
polyTop = [linspace(polyOuter(end, 1),...
            polyOuter(end, 1) - (radiusOuter - radiusInner), geoGenParams.nTop)', polyOuter(end, 2) * ones(geoGenParams.nTop, 1)];

polyInner = m + radiusInner * [cos(phi(end : -1 : 1))', sin(phi(end : -1 : 1))'];
polyLeft = [polyInner(end, 1) * ones(geoGenParams.nTop, 1),...
            linspace(polyInner(end, 2), polyInner(end, 2) - (radiusOuter - radiusInner), geoGenParams.nLeft)'];

geoGenParams.boundaryPolygon = [polyOuter;
                                polyTop(2 : end, :);
                                polyInner(2 : end, :);
                                polyLeft(2 : end, :)];                      

lc = 0.025 * ones(1, size(geoGenParams.boundaryPolygon, 1));

generateGmshGeoFile(geoGenParams.boundaryPolygon, meshName, lc, geoGenParams)
end