function [par] = loadTestCase5(g)
% LOADTESTCASE5 Define the setup for the test case 5.
% 
%   The Gmsh geometry (rectangle) has the following boundary IDs:
%       1: right
%       2: top
%       3: left
%       4: bottom
%   The boundary parameters are stored in different data structures. They
%   are ordered corresponding to the order of the boundary IDs above.
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%
%% domain
par.spatialDimension = 2;
par.meshPrefix = 'meshTestCase5';
% boundaryID: column in outerNormals3D
par.outerNormalVectors3D = [1, 0, 0;
                            0, 1, 0;
                            -1, 0, 0;
                            0, -1, 0]';
                                                                     
%% reflectivity at boundary
% boundaryID: value of reflectivity
% 0: no reflection
% 1: total reflection
par.reflectivity = [0, 0.99, 0, 0.99];

%% scattering kernel
addpath('./testCases/')
par.kernel = henyeyGreenstein(g);

%% external source at boundary
par.externalSource = {@(vx, vy, vz) zeros(size(vx));
                      @(vx, vy, vz) zeros(size(vx));
                      @(vx, vy, vz) vx;
                      @(vx, vy, vz) zeros(size(vx))};
                       
%% absorption, scattering, attenuation coefficient
par.sigma_a = @(x) 0.0 * ones(1, size(x, 2));
par.sigma_s = @(x) 1.0 * ones(1, size(x, 2));
end