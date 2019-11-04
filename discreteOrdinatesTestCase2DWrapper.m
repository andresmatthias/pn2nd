function discreteOrdinatesTestCase2DWrapper(nameTestCase, par, meshNames, loadMesh,...
    maxExactDegree, getBoundaryEdgesIdx2IdTestCase2D_x)
% DISCRETEORDINATESTESTCASE2DWRAPPER Wrapper to run the discrete ordinates
%   method for 2D test cases.
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

buildFolder = pathDeclarationsDiscreteOrdinates(nameTestCase);
createFolder(buildFolder)
meshs = getMesh();

for j = 1 : length(meshs)
    [radiativeEnergyDO, ~, par] = mainDiscreteOrdinates(par, meshs{j}, maxExactDegree(j));
    mesh = meshs{j};
    save(sprintf('%s%s_radEnergyDiscOrdinates_Deg_%d.mat', buildFolder, meshNames{j}, maxExactDegree(j)), 'radiativeEnergyDO', 'mesh', 'par')
end

function meshs = getMesh()
    if loadMesh
        meshs = cell(length(meshNames), 1);
        for k = 1 : length(meshNames)
            meshs{k} = load([buildFolder, meshNames{k}, '.mat']);
        end
    else
        meshs = cell(length(meshNames), 1);
        for k = 1 : length(meshNames)
            meshs{k} = getMeshTestCase2D(buildFolder, meshNames{k}, buildFolder, ...
                getBoundaryEdgesIdx2IdTestCase2D_x);
        end
    end
end
end

function buildFolder = pathDeclarationsDiscreteOrdinates(nameTestCase)
    restoredefaultpath;
    buildFolder = ['./build/', nameTestCase, '/'];
    addpath('./testCases/')
    pathToTestCase = ['./testCases/', nameTestCase, '/'];
    addpath(pathToTestCase);
    addpath('./discreteOrdinates/2D/')
    addpath('./discreteOrdinates/')
    addpath('./tools/')
    addpath('./PN2nd/quadratureRules/')
    addpath('./PN2nd/')
end