function generateTestCase5()
% GENERATETESTCASE5 Generate the second-order formulation of the PN
%   equations for given order and scattering kernel in test case 5.
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
