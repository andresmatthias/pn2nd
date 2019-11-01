function runTestCase3()
% RUNTESTCASE3 Run the test case 3 and visualize the results.
% 
%   User instruction:
%       1) Use generateTestCase3.m first to generate second-order 
%          formulation of PN equations.
%       2) Run runTestCase3.py in Python (with FEniCS) from folder ./FEniCS .
%       3) Then execute this file.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

clear
close all

%% setup
nameTestCase = 'testCase3';
buildFolder = pathDeclarations(nameTestCase);
par = loadTestCase3; 
ModelOrders = [1, 3, 5, 7];
nGridPoints = 1 + [10, 20, 40, 80, 160, 320, 640];
gridDO = linspace(0, 1, 100);
maxExactDegreeDO = 23;


%%
createFolder(buildFolder);
[zGridRef, radEnergyKinetic] = getKineticReferenceSolution(par, nGridPoints);
FEniCSSolution = getFEniCSSolution(buildFolder, nGridPoints, ModelOrders);
radEnergyDO = getDiscreteOrdinatesSolution(buildFolder, par, gridDO, maxExactDegreeDO);
numConvErrors = numConvergenceRateInf(ModelOrders, nGridPoints, FEniCSSolution);
visualize()
writeData2File()


%%
function visualize() 
    figure('Name','Test case 3','NumberTitle','off')
    subplot(1,2,1)
    plot(zGridRef, radEnergyKinetic, 'g+')
    hold on

    plot(FEniCSSolution{1, end}.points, FEniCSSolution{1, end}.radiativeEnergy, 'b')
    plot(FEniCSSolution{2, end}.points, FEniCSSolution{2, end}.radiativeEnergy, 'y')
    plot(FEniCSSolution{3, end}.points, FEniCSSolution{3, end}.radiativeEnergy, 'r')
    plotPwConstantFunction1D(radEnergyDO, gridDO, 'Color', 'k')
    legend('kinetic', 'P1_{2nd}', 'P3_{2nd}', 'P5_{2nd}', 'disc. ord.')
    title('various solutions')

    
    subplot(1,2,2)
    loglog(nGridPoints, [nan, numConvErrors(1, :)])
    title('||P1_{2nd,h} - P1_{2nd, 2h} ||_{\infty} , dx \rightarrow 0')
end

function writeData2File()
    %% kinetic refernce solution
    fileName = [buildFolder, 'kineticReference.csv'];
    headers = {'z', 'radEnergy'};
    matrix = [zGridRef', radEnergyKinetic'];
    write2File(fileName, headers, matrix)
    
    %% PN_2nd
    fileName = [buildFolder, 'PN2nd.csv'];
    headers = {'z', 'N1', 'N3', 'N5'};
    matrix = [zGridRef', FEniCSSolution{1, end}.radiativeEnergy',...
                         FEniCSSolution{2, end}.radiativeEnergy',...
                         FEniCSSolution{3, end}.radiativeEnergy'];
    write2File(fileName, headers, matrix)
    
    %% discrete ordinates
    fileName = [buildFolder, 'discOrd.csv'];
    headers = cell(0,1);
    matrix = zeros(2, (length(gridDO) - 1)*2);
    for k = 1 : length(gridDO) - 1
        matrix(:, (k - 1) * 2 + 1) = [gridDO(k); gridDO(k + 1)];
        matrix(:, (k - 1) * 2 + 2) = [radEnergyDO(k); radEnergyDO(k)];
        headers{end + 1} = sprintf('cellLR%d', k);
        headers{end + 1} = sprintf('radEnergy%d', k);
    end
    write2File(fileName, headers, matrix)
    
    %% numerical convergence P1_2nd
    fileName = [buildFolder, 'numConvErrorsP1_2nd.csv'];
    headers = {'nGridPoints', 'numConvErrorP1_2nd'};
    matrix = [nGridPoints', [numConvErrors(1, :)'; nan]];
    write2File(fileName, headers, matrix)
end
end

function convErrorsInf = numConvergenceRateInf(ModelOrders, nGridPoints, FEniCSSolution)
    convErrorsInf = zeros(length(ModelOrders), length(nGridPoints) - 1);
    for i = 1:length(ModelOrders) 
        for k = 2:length(nGridPoints)
            grid1 = linspace(0, 1, nGridPoints(k - 1));
            grid2 = linspace(0, 1, nGridPoints(k));
            tmp1 = interp1(grid1, FEniCSSolution{i, k - 1}.radiativeEnergy, grid2);

            convErrorsInf(i, k - 1) = max(abs(tmp1 - FEniCSSolution{i, k}.radiativeEnergy));             
        end
    end
    
    for i = 1 : length(ModelOrders)
        fprintf('\n\nnum. convergence (dx -> 0) rate P%d2nd solution \n', ModelOrders(i))
        for k = 2 : length(nGridPoints) - 1
            fprintf('%1.2f\t', log(convErrorsInf(i, k - 1) / convErrorsInf(i, k)) / log(2))
        end
    end
end

function buildFolder = pathDeclarations(nameTestCase)
    restoredefaultpath;
    buildFolder = ['./build/', nameTestCase, '/'];
    pathToTestCase = ['./testCases/', nameTestCase, '/'];
    addpath(pathToTestCase)
    addpath('./testCases/')
    addpath(genpath('./PN2nd/')) 
    addpath('./discreteOrdinates/1D/')
    addpath('./discreteOrdinates/')
    addpath('./tools/')
end

function [zGridRef, radEnergyKinetic] = getKineticReferenceSolution(par, nGridPoints)
zGridRef = linspace(par.meshZMin, par.meshZMax, nGridPoints(end));
radEnergyKinetic = radiativeEnergyKineticSolutionTestCase3(zGridRef);
end

function FEniCSSolution = getFEniCSSolution(buildFolder, nGridPoints, ModelOrders)
    FEniCSSolution = cell(length(ModelOrders), length(nGridPoints));
    for i = 1:length(ModelOrders)
        for j = 1 : length(nGridPoints)
            FEniCSSolution{i, j} = load([buildFolder, 'FEniCS/', sprintf('radiativeEnergy_P%d2nd_n_%d.mat', ModelOrders(i), nGridPoints(j))]);
        end
    end
end

function radEnergyDO = getDiscreteOrdinatesSolution(buildFolder, par, gridDO, maxExactDegreeDO)
    par.grid = gridDO;
    par.nGridPoints = length(par.grid);
    par.maxExactDegree = maxExactDegreeDO;
    radEnergyDO = mainDiscreteOrdinates1D(par);
    save([buildFolder, sprintf('discreteOrdinates_maxDeg_%d.mat', maxExactDegreeDO)], 'radEnergyDO', 'gridDO', 'par');
end

