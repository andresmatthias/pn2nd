function runTestCase2()
% RUNTESTCASE2 Run the test case 2 and visualize the results.
% 
%   User instruction:
%       1) Use generateTestCase2.m first to generate second-order 
%          formulation of PN equations.
%       2) Run runTestCase2.py in Python (with FEniCS) from folder ./FEniCS .
%       3) Then execute this file.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

clear
close all
%% setup
nameTestCase = 'testCase2';
buildFolder = pathDeclarations(nameTestCase);
par = loadTestCase2; 
ModelOrders = 1 : 2 : 21;
nGridPoints = 1 + [10, 20, 40, 80, 160, 320, 640];
gridDO = linspace(0, 1, nGridPoints(end));
gridDOPlot = linspace(0, 1, 100);
maxExactDegreeDO = 23;

loadPNAna = true;

%%
createFolder(buildFolder);
radEnergyPNOrig = getPNOrigSolution(buildFolder, par, nGridPoints, ModelOrders(1 : 5), loadPNAna);
FEniCSSolution = getFEniCSSolution(buildFolder, nGridPoints, ModelOrders);
diffInfPN2FEniCS = getDifferenceInfPNOrig2FEniCSSolution(radEnergyPNOrig, FEniCSSolution, nGridPoints, ModelOrders(1 : 5));
radEnergyDO = getDiscreteOrdinatesSolution(buildFolder, par, gridDO, maxExactDegreeDO);
radEnergyDOPlot = getDiscreteOrdinatesSolution(buildFolder, par, gridDOPlot, maxExactDegreeDO);
numConvErrors = numConvergenceRateInf(ModelOrders, nGridPoints, FEniCSSolution); 
diffInfRelDO2PN2nd = getDiffInfRelDO2PN2nd(ModelOrders, FEniCSSolution, radEnergyDO);
visualize()
writeData2File()

%%
function visualize() 
    figure('Name','Test case 2','NumberTitle','off')
    subplot(1, 3, 1)
    plot(FEniCSSolution{1, end}.points, FEniCSSolution{1, end}.radiativeEnergy, 'b')
    hold on
    plot(FEniCSSolution{2, end}.points, FEniCSSolution{2, end}.radiativeEnergy, 'g')
    plot(FEniCSSolution{end, end}.points, FEniCSSolution{end, end}.radiativeEnergy, 'r')
    plotPwConstantFunction1D(radEnergyDO, gridDO, 'Color', 'k')
    legend('P1_{2nd}', 'P3_{2nd}', 'P21_{2nd}', 'DO')
    title('various solutions')

    subplot(1, 3, 2)
    for k = 3:length(ModelOrders)
       fprintf('\nloglog conv Error DO vs PN2nd %1.2f', log((diffInfRelDO2PN2nd(k-2) - diffInfRelDO2PN2nd(k-1)) / (diffInfRelDO2PN2nd(k-1) - diffInfRelDO2PN2nd(k))) / log(2));
    end
    semilogy(ModelOrders, diffInfRelDO2PN2nd)
    title('||PN2nd - DO||_{\infty} / ||DO||_{\infty}, n\rightarrow \infty')

    subplot(1, 3, 3)
    for k = 1 : size(radEnergyPNOrig, 1)
        loglog(nGridPoints, diffInfPN2FEniCS(k, :))
        hold on
    end 
    title('||PN-PN2nd||_{\infty}, dx \rightarrow 0')
    for j = 1 : size(radEnergyPNOrig, 1)
        fprintf('\n\nnum. convergence (dx -> 0) rate P%d2nd solution \n', ModelOrders(j))
        for k = 2:length(nGridPoints)
            fprintf('%1.2f\t', log(diffInfPN2FEniCS(j, k-1) / diffInfPN2FEniCS(j, k)) / log(2));
        end
    end
end


function writeData2File()
    zGridFine = linspace(0, 1, nGridPoints(end));
    %% PN2nd
    fileName = [buildFolder, 'PN2nd.csv'];
    headers = cell(length(ModelOrders) + 1, 1);
    headers{1} = 'z';
    headers(2:end) = compose('N%d', ModelOrders);
    matrix = zGridFine';
    for k = 1 : length(ModelOrders)
       matrix(:, end + 1) = FEniCSSolution{k, end}.radiativeEnergy'; 
    end
    write2File(fileName, headers, matrix)
    
    %% discrete ordinates
    fileName = [buildFolder, 'discOrd.csv'];
    headers = cell(0,1);
    matrix = zeros(2, (length(gridDOPlot) - 1)*2);
    for k = 1 : length(gridDOPlot) - 1
        matrix(:, (k - 1) * 2 + 1) = [gridDOPlot(k); gridDOPlot(k + 1)];
        matrix(:, (k - 1) * 2 + 2) = [radEnergyDOPlot(k); radEnergyDOPlot(k)];
        headers{end + 1} = sprintf('cellLR%d', k);
        headers{end + 1} = sprintf('radEnergy%d', k);
    end
    write2File(fileName, headers, matrix)
    
    %% numerical convergence PN2nd
    fileName = [buildFolder, 'numConvErrorsPN2nd.csv'];
    headers = cell(length(ModelOrders) + 1, 1);
    headers{1} = 'nGridPoints';
    headers(2:end) = compose('diff2CoarseN%d', ModelOrders);
    matrix = nGridPoints';
    for k = 1 : length(ModelOrders)
        matrix(:, end + 1) = [nan; numConvErrors(k, :)'];
    end
    write2File(fileName, headers, matrix)
    
    %% difference PN2nd 2 DO
    fileName = [buildFolder, 'diffInfRelDO2PN2nd.csv'];
    headers = {'ModelOrders', 'diffInfRelDO2PN2nd'};
    matrix = [ModelOrders', diffInfRelDO2PN2nd'];
    write2File(fileName, headers, matrix)
    
    %% difference (inf) PN original 2 PN2nd (FEniCS) 
    fileName = [buildFolder, 'diffInfPNOrig2PN2nd.csv'];
    n = size(diffInfPN2FEniCS, 1);
    headers = cell(n + 1, 1);
    headers{1} = 'nGridPoints';
    headers(2:end) = compose('N%d', 1 : 2 : 2 * n);
    matrix = [nGridPoints',diffInfPN2FEniCS'];
    write2File(fileName, headers, matrix)
end
end

function diffDO2PN2nd = getDiffInfRelDO2PN2nd(ModelOrders, FEniCSSolution, radEnergyDO)
    diffDO2PN2nd = zeros(size(ModelOrders));
    for k = 1 : length(ModelOrders)
        e1 = max(abs( FEniCSSolution{k, end}.radiativeEnergy(1 : end-1) - radEnergyDO));
        e2 = max(abs(FEniCSSolution{k, end}.radiativeEnergy(2 : end) - radEnergyDO));
       diffDO2PN2nd(k) = max(e1, e2) / max(abs(radEnergyDO)); 
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


function radEnergyPNOrig = getPNOrigSolution(buildFolder,  par, nGridPoints, ModelOrders, loadPNAna)
    options.RelTol = 1e-10;
    options.AbsTol = 1e-10;
    if loadPNAna
        load([buildFolder, 'radEnergyPNOrig.mat'], 'radEnergyPNOrig');
    else
        radEnergyPNOrig = cell(length(ModelOrders), length(nGridPoints));
        for i = 1 : length(ModelOrders)
            for j = 1 : length(nGridPoints)
                fprintf('compute solution of original PN equations: %d / %d, %d / %d\n', i, length(ModelOrders), j, length(nGridPoints))
                radEnergyPNOrig{i, j} = radiativeEnergyPNOrigIsotropicKernel1D(par, ModelOrders(i), nGridPoints(j), options); 
            end
        end
    save([buildFolder, 'radEnergyPNOrig.mat'], 'radEnergyPNOrig')
    end
end

function FEniCSSolution = getFEniCSSolution(buildFolder, nGridPoints, ModelOrders)
    FEniCSSolution = cell(length(ModelOrders), length(nGridPoints));
    for i = 1:length(ModelOrders)
        for j = 1 : length(nGridPoints)
            FEniCSSolution{i, j} = load([buildFolder, 'FEniCS/', sprintf('radiativeEnergy_P%d2nd_n_%d.mat', ModelOrders(i), nGridPoints(j))]);
        end
    end
end

function [diffInfPNOrig2FEniCS] = getDifferenceInfPNOrig2FEniCSSolution(radEnergyPNOrig, FEniCSSolution, nGridPoints, ModelOrders)
    diffInfPNOrig2FEniCS = zeros(length(ModelOrders), length(nGridPoints));
    for i = 1:length(ModelOrders)
        for j = 1 : length(nGridPoints)
           diffInfPNOrig2FEniCS(i, j) = max(abs(radEnergyPNOrig{i, j} - FEniCSSolution{i, j}.radiativeEnergy ));
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

