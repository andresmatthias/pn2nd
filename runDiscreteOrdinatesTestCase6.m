function runDiscreteOrdinatesTestCase6()
% RUNDISCRETEORDINATESTESTCASE6 Compute the discrete ordinates 
%   approximation for the solution of the kinetic problem of test case 6.
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

clear
close all
%% setup
nameTestCase = 'testCase6';
pathToTestCase = ['./testCases/', nameTestCase, '/'];
addpath(pathToTestCase);
par = loadTestCase6; 
meshNames = {'meshTestCase6_0', 'meshTestCase6_0'};
loadMesh = false;
aux = load(sprintf('./build/%s/meshTestCase6Normals.mat',nameTestCase));
getBoundaryEdgesIdx2IdTestCase2D_x =...
    @(TR, points, edgeList) getBoundaryEdgesIdx2IdTestCase6(TR, points, edgeList, aux.geoGenParams);
maxExactDegree = [15, 23];

%%
discreteOrdinatesTestCase2DWrapper(nameTestCase, par, meshNames, loadMesh,...
    maxExactDegree, getBoundaryEdgesIdx2IdTestCase2D_x)
end