function generateTestCase2()
% GENERATETESTCASE2 Generate the second-order formulation of the PN
%   equations for given order and scattering kernel in test case 2.
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
