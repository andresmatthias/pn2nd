function generateTestCase4()
% GENERATETESTCASE4 Generate the second-order formulation of the PN
%   equations for given order and scattering kernel in test case 4.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

clear
close all
%% setup
nameTestCase = 'testCase4';
loadSystemDomain = false;
loadSystemBoundary = false;
ModelOrders = [1, 3, 5, 7];
pathToTestCase = ['./testCases/', nameTestCase, '/'];
addpath(pathToTestCase)
par = loadTestCase4;
%%
generatePN2nd(nameTestCase, ModelOrders, loadSystemDomain, loadSystemBoundary, par)
end
