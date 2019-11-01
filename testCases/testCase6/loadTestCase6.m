function [par] = loadTestCase6()
% LOADTESTCASE6 Define the setup for the test case 6.
% 
%   The Gmsh geometry (rectangle) has the following boundary IDs:
%       1: right
%       2: top
%       3: left
%       4: bottom
%   The boundary parameters are stored in different data structures. They
%   are ordered corresponding to the order of the boundary IDs above.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%
%% domain
par.spatialDimension = 2;
par.meshPrefix = 'meshTestCase6';
% boundaryID: column in outerNormals3D
aux = load('./build/testCase6/meshTestCase6Normals.mat');
par.outerNormalVectors3D = aux.outerNormalVectors3D;
% the following parameters are defined in generateTestCase6GeoFile.m
par.geoGenParams = aux.geoGenParams;
nTop = aux.geoGenParams.nTop;
nInnerOuter = aux.geoGenParams.nInnerOuter;
nLeft = aux.geoGenParams.nLeft;
                                                                     
%% reflectivity at boundary
% boundaryID: value of reflectivity
% 0: no reflection
% 1: total reflection
par.reflectivity = [0.5 * ones(1, nInnerOuter - 1), zeros(1, nTop - 1),...
                    0.5 * ones(1, nInnerOuter - 1), zeros(1, nLeft - 1)];

%% scattering kernel
addpath('./testCases/')
par.kernel = isotropicKernel;

%% external source at boundary
par.externalSource = cell(2 * (nInnerOuter - 1) + (nTop - 1) + (nLeft - 1), 1);
%% outer curve
for k = 1 : nInnerOuter - 1
    par.externalSource{k} = @(vx, vy, vz) zeros(size(vx));
end
%% top
for k = 1 : nTop - 1
    par.externalSource{k + nInnerOuter - 1} = @(vx, vy, vz) zeros(size(vx));
end
%% inner curve
for k = 1 : nInnerOuter - 1
    par.externalSource{k + nInnerOuter + nTop - 2} = @(vx, vy, vz) zeros(size(vx));
end
%% left
for k = 1 : nLeft - 1
    par.externalSource{k + nInnerOuter *2 + nTop - 3} = @(vx, vy, vz) 1.0 * ones(size(vx));
end

%% absorption, scattering, attenuation coefficient
par.sigma_a = @(x) 0.0 * ones(1, size(x, 2));
par.sigma_s = @(x) 0.1 * ones(1, size(x, 2));
end