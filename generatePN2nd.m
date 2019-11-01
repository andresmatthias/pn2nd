function generatePN2nd(nameTestCase, ModelOrders, loadSystemDomain, ...
    loadSystemBoundary, par, varargin)
% GENERATEPN2ND Generate the second-order formulation of the PN-equations.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

if nargin == 5
    [buildFolder, pathToPN2ndModels] = pathDeclarations(nameTestCase);
else
    buildFolder = varargin{1};
    pathToPN2ndModels = varargin{2};
end
createFolder(buildFolder);
createFolder(pathToPN2ndModels);

for N = ModelOrders
    %% 
    fprintf('generate second-order formulation of PN equations on domain for specific kernel: %d / %d\n', N, ModelOrders(end))
    if loadSystemDomain
        systemMatricesDomain = load(sprintf('%ssystemMatricesDomain_d_%d_k_%s_N_%d.mat', pathToPN2ndModels,...
            par.spatialDimension, par.kernel.name, N));
    else
        systemMatricesDomain = generatePN2ndWeakSystemOnDomain(N, par.kernel, par.spatialDimension);
        save(sprintf('%ssystemMatricesDomain_d_%d_k_%s_N_%d.mat', pathToPN2ndModels,...
            par.spatialDimension, par.kernel.name, N), '-struct', 'systemMatricesDomain');
    end

    systemMatricesDomain2Python(N, par.spatialDimension, par.kernel.name,...
        pathToPN2ndModels, systemMatricesDomain)

    %% 
    fprintf('generate second-order formulation of PN equations on boundary for specific outer normals\n')
    if loadSystemBoundary 
        systemBoundary = load(sprintf('%ssystemBoundary_%s_N_%d.mat', buildFolder,...
            nameTestCase, N));
    else
        systemBoundary = generatePN2ndBoundarySystem(par.outerNormalVectors3D, N, par.spatialDimension);
        save(sprintf('%ssystemBoundary_%s_N_%d.mat', buildFolder,...
            nameTestCase, N), '-struct', 'systemBoundary');
    end
end
end

function [buildFolder, pathToPN2ndModels] = pathDeclarations(nameTestCase)
    restoredefaultpath
    buildFolder = ['./build/', nameTestCase, '/'];
    pathToPN2ndModels = './build/PN2ndModels/';
    pathToTestCase = ['./testCases/', nameTestCase, '/'];
    addpath(pathToTestCase)
    addpath('./testCases/')
    addpath('./tools/')
    addpath(genpath('./PN2nd/'))
end