function runDiscreteOrdinatesTestCase5()
% RUNDISCRETEORDINATESTESTCASE5 Compute the discrete ordinates 
%   approximation for the solution of the kinetic problem of test case 5.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
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