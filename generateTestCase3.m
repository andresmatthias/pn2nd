function generateTestCase3()
clear
close all
%% setup
nameTestCase = 'testCase3';
loadSystemDomain = false;
loadSystemBoundary = false;
ModelOrders = 1 : 2 : 7;
pathToTestCase = ['./testCases/', nameTestCase, '/'];
addpath(pathToTestCase)
par = loadTestCase3;
%%
generatePN2nd(nameTestCase, ModelOrders, loadSystemDomain, loadSystemBoundary, par)
end
