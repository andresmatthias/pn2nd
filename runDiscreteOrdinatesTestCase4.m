function runDiscreteOrdinatesTestCase4()
% RUNDISCRETEORDINATESTESTCASE4 Compute the discrete ordinates 
%   approximation for the solution of the kinetic problem of test case 4.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

clear
close all
%% setup
nameTestCase = 'testCase4';
pathToTestCase = ['./testCases/', nameTestCase, '/'];
addpath(pathToTestCase);
par = loadTestCase4; 
meshNames = {'meshTestCase4_0', 'meshTestCase4_0'};
loadMesh = false;
getBoundaryEdgesIdx2IdTestCase2D_x =...
    @(TR, points, edgeList) getBoundaryEdgesIdx2IdTestCase4(TR, points, edgeList);
maxExactDegree = [15, 23];

%%
discreteOrdinatesTestCase2DWrapper(nameTestCase, par, meshNames, loadMesh,...
    maxExactDegree, getBoundaryEdgesIdx2IdTestCase2D_x)
end