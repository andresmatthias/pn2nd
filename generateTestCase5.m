function generateTestCase5()
% GENERATETESTCASE5 Generate the second-order formulation of the PN
%   equations for given order and scattering kernel in test case 5.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%
clear
close all
%% setup
g = [0.0, 0.5];
for gi = g
    nameTestCase = 'testCase5';
    loadSystemDomain = false;
    loadSystemBoundary = false;
    ModelOrders = [1, 3, 5, 7];
    pathToTestCase = ['./testCases/', nameTestCase, '/'];
    addpath(pathToTestCase)
    par = loadTestCase5(gi);
    %%
    generatePN2nd(nameTestCase, ModelOrders, loadSystemDomain, loadSystemBoundary, par)
end
end
