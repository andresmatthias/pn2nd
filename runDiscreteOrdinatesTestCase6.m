function runDiscreteOrdinatesTestCase6()
% RUNDISCRETEORDINATESTESTCASE6 Compute the discrete ordinates 
%   approximation for the solution of the kinetic problem of test case 6.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
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