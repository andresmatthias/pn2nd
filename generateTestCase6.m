function generateTestCase6()
% GENERATETESTCASE6 Generate the second-order formulation of the PN
%   equations for given order and scattering kernel in test case 6.
% 
%   User instruction:
%       1) Execute generateTestCase6GeoFile.m in folder
%           ./testCases/testCase6 to create the .geo file for this
%           customized domain.
%       2) Execute genMeshTestCase6.sh in folder ./testCases/testCase6 in
%           the command line (knowing dolfin and gmsh).
%       3) Then execute this file.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%
clear
close all
%% setup
nameTestCase = 'testCase6';
loadSystemDomain = false;
loadSystemBoundary = false;
ModelOrders = [1, 3, 5, 7];
pathToTestCase = ['./testCases/', nameTestCase, '/'];
addpath(pathToTestCase)
par = loadTestCase6;
%%
generatePN2nd(nameTestCase, ModelOrders, loadSystemDomain, loadSystemBoundary, par)

end
