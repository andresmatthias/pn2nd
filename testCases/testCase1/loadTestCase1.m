function [par] = loadTestCase1()
% LOADTESTCASE1 Define the setup for the test case 1.
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
par.outerNormalVectors3D = [0, 0, -1;
                            0, 0, 1]';

%% reflectivity at boundary
% 0: no reflection
% 1: total reflection
par.reflectivity = [0, 0];

%% scattering kernel
addpath('./testCases/')
par.kernel = isotropicKernel;

%% external source at boundary
par.externalSource = {@(vx, vy, vz) 1 / 4 / pi * ones(size(vx));
                      @(vx, vy, vz) zeros(size(vx))};

%% absorption, scattering, attenuation coefficient
par.sigma_a = @(z) ones(1, size(z, 2));
par.sigma_s = @(z) zeros(1, size(z, 2));
end