function runTestCase1()
% RUNTESTCASE1 Run the test case 1 and visualize the results.
% 
%   User instruction:
%       1) Use generateTestCase1.m first to generate second-order 
%          formulation of PN equations.
%       2) Run runTestCase1.py in Python (with FEniCS) in folder ./FEniCS .
%       3) Then execute this file.
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
nameTestCase = 'testCase1';
buildFolder = pathDeclarations(nameTestCase);
par = loadTestCase1; 
ModelOrders = 1 : 2 : 21;
nGridPoints = 1 + [10, 20, 40, 80, 160, 320, 640];
gridDO = linspace(0, 1, 100);
maxExactDegreeDO = 23;

loadPNAna = true;

%%
createFolder(buildFolder)
[zGridRef, radEnergyKinetic] = getKineticReferenceSolution(par, nGridPoints);
radEnergyPNAnalytic = getPNAnalyticSolution(buildFolder, par, nGridPoints, ModelOrders, loadPNAna);
[diffInfPNAnalytic2FEniCS, ~] = getDifferenceInfPNAnalytic2FEniCSSolution(buildFolder, radEnergyPNAnalytic, nGridPoints, ModelOrders);
[diffInfRelKin2PNAna] = getRelativeDifferenceKin2PNAna(radEnergyKinetic, radEnergyPNAnalytic, ModelOrders);
radEnergyDO = getDiscreteOrdinatesSolution(buildFolder, par, gridDO, maxExactDegreeDO);
visualize()
writeData2File()

%%
function writeData2File()
    %% relative difference (inf) kinetic solution 2 analytical PN solution
    fileName = [buildFolder, 'diffInfRelKin2PNAna.csv'];
    headers = {'ModelOrders', 'diffInfRelKin2PNAna'};
    matrix = [ModelOrders', diffInfRelKin2PNAna'];
    write2File(fileName, headers, matrix)
    
    %% difference (inf) PN analytic 2 FEniCS 
    fileName = [buildFolder, 'diffInfPNAnalytic2FEniCS.csv'];
    headers = cell(length(ModelOrders) + 1, 1);
    headers{1} = 'nGridPoints';
    headers(2:end) = compose('N%d', ModelOrders);
    matrix = [nGridPoints',diffInfPNAnalytic2FEniCS'];
    write2File(fileName, headers, matrix)
    
    %% kinetic refernce solution
    fileName = [buildFolder, 'kineticReference.csv'];
    headers = {'z', 'radEnergy'};
    matrix = [zGridRef', radEnergyKinetic];
    write2File(fileName, headers, matrix)
    
    %% P1 analytic
    fileName = [buildFolder, 'P1Ana.csv'];
    headers = {'z', 'radEnergy'};
    matrix = [zGridRef', radEnergyPNAnalytic{1, end}'];
    write2File(fileName, headers, matrix)
    
    %% PMax analytic
    fileName = [buildFolder, sprintf('P%dAna.csv',ModelOrders(end))];
    headers = {'z', 'radEnergy'};
    matrix = [zGridRef', radEnergyPNAnalytic{end, end}'];
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
end

%%
function visualize() 
    figure('Name','Test case 1','NumberTitle','off')
    subplot(1, 3, 1)
    plot(zGridRef, radEnergyKinetic, 'k-')
    hold on
    plot(zGridRef, radEnergyPNAnalytic{1, end})
    plot(zGridRef, radEnergyPNAnalytic{end, end})
    plotPwConstantFunction1D(radEnergyDO, gridDO, 'Color', 'k')
    legend('Kin', 'Ana P1', 'Ana P21', 'DO')
    title('various solutions')

    subplot(1, 3, 2)
    semilogy(ModelOrders, diffInfRelKin2PNAna)
    title('||PN - kin||_{\infty} / ||kin||_{\infty}, N\rightarrow \infty')

    subplot(1, 3, 3)
    for i = 1 : length(ModelOrders)
        loglog(nGridPoints, diffInfPNAnalytic2FEniCS(i, :))
        hold on
    end
    title('||PN - PN2nd||_{\infty}, dx \rightarrow 0')
    for i = 1:length(ModelOrders)
        fprintf('\n\nnum. convergence (dx -> 0) rate SP%d solution \n', ModelOrders(i))
        for k = 2:length(nGridPoints)
            fprintf('%1.2f\t', log(diffInfPNAnalytic2FEniCS(i, k-1) / diffInfPNAnalytic2FEniCS(i, k)) / log(2));
        end
    end
end
end

%%
function [diffInfRelKin2PNAna] = getRelativeDifferenceKin2PNAna(...
    radEnergyKinetic, radEnergyPNAnalytic,ModelOrders)
% Compute the relative difference between analytic solution of PN 
% equations and solution of kinetic equation. 
    diffInfRelKin2PNAna = zeros(size(ModelOrders));
    for k = 1 : length(ModelOrders)
        diffInfRelKin2PNAna(k) = max(abs(radEnergyPNAnalytic{k, end}') - radEnergyKinetic)...
                                    / max(abs(radEnergyKinetic)); 
    end
end

%%
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

%%
function [zGridRef, radEnergyKinetic] = getKineticReferenceSolution(par, nGridPoints)
% Compute the solution of the kinetic equation.
    zGridRef = linspace(par.meshZMin, par.meshZMax, nGridPoints(end));
    radEnergyKinetic = radiativeEnergyKineticSolutionTestCase1(par.externalSource{1},...
        par.externalSource{2}, zGridRef, par.sigma_a(1));
end

%%
function radEnergyPNAnalytic = getPNAnalyticSolution(buildFolder, par, nGridPoints, ModelOrders, loadPNAna)
% Compute the solution of the original (1st order) PN equations.
    if loadPNAna
        load([buildFolder, 'radEnergyPNAnalytic.mat'], 'radEnergyPNAnalytic');
    else
        radEnergyPNAnalytic = cell(length(ModelOrders), length(nGridPoints));
    for i = 1 : length(ModelOrders)
        for j = 1 : length(nGridPoints)
            fprintf('compute analytic PN solution: %d / %d, %d / %d\n', i, length(ModelOrders), j, length(nGridPoints))
            radEnergyPNAnalytic{i, j} = radiativeEnergyPNOrigIsotropicKernel1D(par, ModelOrders(i), nGridPoints(j), 'homogeneous'); 
        end
    end
    save([buildFolder, 'radEnergyPNAnalytic.mat'], 'radEnergyPNAnalytic')
    end
end

%%
function [diffInfPNAnalytic2FEniCS, FEniCSSolution] = getDifferenceInfPNAnalytic2FEniCSSolution(buildFolder, radEnergyPNAnalytic, nGridPoints, ModelOrders)
    diffInfPNAnalytic2FEniCS = zeros(length(ModelOrders), length(nGridPoints));
% Compute the inf-difference between the solutions of the original PN
% equations (1st order) and the solution of the second order formulation, 
% compute with the help of the Python toolbox FEniCS.
    FEniCSSolution = cell(length(ModelOrders), length(nGridPoints));
    for i = 1:length(ModelOrders)
        for j = 1 : length(nGridPoints)
            FEniCSSolution{i, j} = load([buildFolder, 'FEniCS/', sprintf('radiativeEnergy_P%d2nd_n_%d.mat', ModelOrders(i), nGridPoints(j))]);
            diffInfPNAnalytic2FEniCS(i, j) = max(abs(radEnergyPNAnalytic{i, j} - FEniCSSolution{i, j}.radiativeEnergy ));
        end
    end
end

%%
function radEnergyDO = getDiscreteOrdinatesSolution(buildFolder, par, gridDO, maxExactDegreeDO)
% Approximate the solution of the kinetic equation with the discrete
% ordinates method.
    par.grid = gridDO;
    par.nGridPoints = length(par.grid);
    par.maxExactDegree = maxExactDegreeDO;
    radEnergyDO = mainDiscreteOrdinates1D(par);
    save([buildFolder, sprintf('discreteOrdinates_maxDeg_%d.mat', maxExactDegreeDO)], 'radEnergyDO', 'gridDO', 'par');
end

