function runDiscreteOrdinatesTestCase5()
% RUNDISCRETEORDINATESTESTCASE5 Compute the discrete ordinates 
%   approximation for the solution of the kinetic problem of test case 5.
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
g = 0.5;
nameTestCase = 'testCase5';
pathToTestCase = ['./testCases/', nameTestCase, '/'];
addpath(pathToTestCase);
par = loadTestCase5(g); 
meshNames = {'meshTestCase5_0', 'meshTestCase5_0'};
loadMesh = false;
getBoundaryEdgesIdx2IdTestCase2D_x =...
    @(TR, points, edgeList) getBoundaryEdgesIdx2IdTestCase5(TR, points, edgeList);
maxExactDegree = [15, 23];

%%
discreteOrdinatesTestCase2DWrapper(nameTestCase, par, meshNames, loadMesh,...
    maxExactDegree, getBoundaryEdgesIdx2IdTestCase2D_x)
end