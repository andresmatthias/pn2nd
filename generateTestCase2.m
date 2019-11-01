function generateTestCase2()
% GENERATETESTCASE2 Generate the second-order formulation of the PN
%   equations for given order and scattering kernel in test case 2.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

clear
close all
%% setup
nameTestCase = 'testCase2';
loadSystemDomain = false;
loadSystemBoundary = false;
ModelOrders = 1 : 2 : 21;
pathToTestCase = ['./testCases/', nameTestCase, '/'];
addpath(pathToTestCase)
par = loadTestCase2;
%%
generatePN2nd(nameTestCase, ModelOrders, loadSystemDomain, loadSystemBoundary, par)
end
