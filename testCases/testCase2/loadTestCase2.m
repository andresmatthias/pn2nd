function [par] = loadTestCase2()
% LOADTESTCASE2 Define the setup for the test case 2.
%
% The first entry in the structure of each boundary parameter corresponds 
% to zMin (left 1D boundary), the second entry corresponds to zMax (right
% 1D boundary).
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

%% domain
par.spatialDimension = 1;
par.meshZMin = 0;
par.meshZMax = 1;
% boundaryID: column in outerNormals3D
par.outerNormalVectors3D = [0, 0, -1;
                            0, 0, 1]';

%% reflectivity at boundary
% boundaryID: value of reflectivity
% 0: no reflection
% 1: total reflection
par.reflectivity = [0.5, 0.5];

%% scattering kernel
addpath('./testCases/')
par.kernel = isotropicKernel;

%% external source at boundary
par.externalSource = {@(vx, vy, vz) vz.^2;
                      @(vx, vy, vz) zeros(size(vx))};

%% absorption, scattering, attenuation coefficient
par.sigma_a = @(z) (2 + sin(2 * pi * z)) / 10;
par.sigma_s = @(z) (3 - z.^2) / 10;
end