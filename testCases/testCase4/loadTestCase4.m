function [par] = loadTestCase4()
% LOADTESTCASE4 Define the setup for the test case 4.
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
par.meshPrefix = 'meshTestCase4';
% boundaryID: column in outerNormals3D
par.outerNormalVectors3D = [1, 0, 0;
                            0, 1, 0;
                            -1, 0, 0;
                            0, -1, 0]';
                                                                     
%% reflectivity at boundary
% boundaryID: value of reflectivity
% 0: no reflection
% 1: total reflection
par.reflectivity = [0, 0.5, 0, 0.5];

%% scattering kernel
addpath('./testCases/')
par.kernel = isotropicKernel;

%% external source at boundary
par.externalSource = {@(vx, vy, vz) zeros(size(vx));
                      @(vx, vy, vz) zeros(size(vx));
                      @(vx, vy, vz) 1 / 4 / pi * ones(size(vx));
                      @(vx, vy, vz) zeros(size(vx))};
                       
%% absorption, scattering, attenuation coefficient

par.sigma_a = @(x) 100 .* f1(x) .* f2(x) .* f3(x);
par.sigma_s = @(x) 0.01 * ones(1, size(x, 2));
end

function out = f1(x)
    idx = (x(1, :) <= 0.6);
    out = ones(1, size(x, 2));
    out(idx) = exp(-(x(1, idx) - 0.6).^2 / 0.01);
end

function out = f2(x)
    idx = (x(1, :) >= 0.7);
    out = ones(1, size(x, 2));
    out(idx) = exp(-(x(1, idx) - 0.7).^2 / 0.01);
end

function out = f3(x)
    idx = (x(2, :) >= 0.4);
    out = ones(1, size(x, 2));
    out(idx) = exp(-(x(2, idx) - 0.4).^2 / 0.01);
end

