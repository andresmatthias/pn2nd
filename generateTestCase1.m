function generateTestCase1()
% GENERATETESTCASE1 Generate the second-order formulation of the PN
%   equations for given order and scattering kernel in test case 1.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

clear
close all
%% setup
nameTestCase = 'testCase1';
loadSystemDomain = false;
loadSystemBoundary = false;
ModelOrders = 1 : 2 : 21;
pathToTestCase = ['./testCases/', nameTestCase, '/'];
addpath(pathToTestCase)
par = loadTestCase1;
%%
generatePN2nd(nameTestCase, ModelOrders, loadSystemDomain, loadSystemBoundary, par)
end
