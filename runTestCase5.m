function runTestCase5()
% RUNTESTCASE5 Run the test case 5 and visualize the results.
% 
%   User instruction:
%       1) Use generateTestCase5.m first to generate second-order 
%          formulation of PN equations.
%       2) Execute genMeshTestCase5.sh in folder ./testCases/testCase5 in
%           the command line (knowing dolfin and gmsh) 
%       3) Run runTestCase5.py in Python (with FEniCS) in folder ./FEniCS .
%       4) Run runDiscreteOrdinatesTestCase5.m
%       5) Then execute this file.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

clear
close all
%% setup
nameTestCaseMain = 'testCase5';
nameTestCase = 'testCase5_g0,500';
nameTestCaseIso = 'testCase5_g0,000';
buildFolderMain = pathDeclarations(nameTestCaseMain);
buildFolder = pathDeclarations(nameTestCase);
buildFolderIso = pathDeclarations(nameTestCaseIso);
meshName = 'meshTestCase5_0';
meshNameFine = 'meshTestCase5_1';
maxExactDegree = [15, 23];

%% discrete ordinates
solutionDOCoarse = load(sprintf('%s%s_radEnergyDiscOrdinates_Deg_%d.mat', buildFolderMain, meshName, maxExactDegree(1)));
solutionDO = load(sprintf('%s%s_radEnergyDiscOrdinates_Deg_%d.mat', buildFolderMain, meshName, maxExactDegree(2)));

%% PN2nd
ax_min = min(solutionDO.radiativeEnergyDO);
ax_max = max(solutionDO.radiativeEnergyDO);
PN_2nd = {};
PN_2ndIso = {};
PN_2ndFine = {};
L2DifferenceRelFineMesh = zeros(4, 1);
L2DifferenceRelDO = zeros(4, 1);
L2DifferenceRelDODOCoarse = L2NormElementFunction(solutionDOCoarse.radiativeEnergyDO - solutionDO.radiativeEnergyDO, solutionDO.mesh.points, solutionDO.mesh.connectivityList) ...
    / L2NormElementFunction(solutionDO.radiativeEnergyDO, solutionDO.mesh.points, solutionDO.mesh.connectivityList);
fprintf('relative L2 difference DO: %1.2e\n', L2DifferenceRelDODOCoarse)

for k = [1, 3, 5, 7]
    PN_2nd{end + 1} = load(sprintf('%sFEniCS/radiativeEnergy_P%d2nd_%s.mat', buildFolder, k, meshName));
    PN_2ndFine{end + 1} = load(sprintf('%sFEniCS/radiativeEnergy_P%d2nd_%s.mat', buildFolder, k, meshNameFine));
    PN_2ndIso{end + 1} = load(sprintf('%sFEniCS/radiativeEnergy_P%d2nd_%s.mat', buildFolderIso, k, meshName));
    [dPN2ndRefInfRel, L2DifferenceRelFineMesh((k+1)/2)] = relativeDistancePN2ndRefFine(PN_2ndFine{end}, PN_2nd{end});
    L2DifferenceRelDO((k+1)/2) = relativeDistancePN2ndDO(PN_2nd{end}, solutionDO);
    fprintf('Relative difference (inf) P%d: ||coarse - fine||_inf / ||fine||_inf: \t%f\n', k, dPN2ndRefInfRel);
    fprintf('Relative difference (L2) P%d: ||coarse - fine||_L2 / ||fine||_L2: \t%f\n', k, L2DifferenceRelFineMesh((k+1)/2));
    fprintf('Relative difference (L2) ||P%d - DO||_L2 / ||DO||_L2: \t\t\t%f\n\n', k, L2DifferenceRelDO((k+1)/2));

    ax_min = min(ax_min, min(PN_2nd{end}.radiativeEnergy));
    ax_max = max(ax_max, max(PN_2nd{end}.radiativeEnergy));
end


%% evaluation along line
mesh.points = PN_2nd{1}.points';
mesh.connectivityList = double(PN_2nd{1}.connectivityList' + 1);
centerLine = [linspace(0, 1, 100); 0.5 * ones(1, 100)];
PN_2ndAlongLine = {};
PN_2ndIsoAlongLine = {};
for k = [1, 3, 5, 7]
    PN_2ndAlongLine{end + 1} = interpMeshFunction(PN_2nd{(k+1)/2}.radiativeEnergy, mesh, centerLine);
    PN_2ndIsoAlongLine{end + 1} = interpMeshFunction(PN_2ndIso{(k+1)/2}.radiativeEnergy, mesh, centerLine);
end
DOAlongLine = evalElementFunction(solutionDO.radiativeEnergyDO, solutionDO.mesh, centerLine); 

visualize()
writeData2File()


function visualize()
    cmap = viridis(256);

    mesh.points = PN_2nd{1}.points';
    mesh.connectivityList = double(PN_2nd{1}.connectivityList' + 1);
    meshFine.points = PN_2ndFine{1}.points';
    meshFine.connectivityList = double(PN_2ndFine{1}.connectivityList' + 1);

    %DO
    figure
    colormap(cmap);
    plotPWConstantFunction2D(mesh, repmat(solutionDO.radiativeEnergyDO, 1, 3), 0)
    colorbar; shading interp; view(2)
    caxis([ax_min, ax_max])
    title('discrete Ordinates')
    fprintf('DO max / min: \t\t %f \t %f\n', max(solutionDO.radiativeEnergyDO), min(solutionDO.radiativeEnergyDO))

    figure
    colormap(cmap);
    plotPWConstantFunction2D(solutionDO.mesh, repmat(solutionDOCoarse.radiativeEnergyDO - solutionDO.radiativeEnergyDO, 1, 3), 0)
    colorbar; shading interp; view(2)
    title('DO coarse - DO')
    fprintf('diff DO max / min: \t %f \t %f\n', max(solutionDOCoarse.radiativeEnergyDO - solutionDO.radiativeEnergyDO),...
                                                  min(solutionDOCoarse.radiativeEnergyDO - solutionDO.radiativeEnergyDO))

    for k = [1, 3, 5, 7]
        figure
        colormap(cmap);
        trisurf(mesh.connectivityList', mesh.points(1, :), mesh.points(2, :), PN_2nd{(k+1)/2}.radiativeEnergy);
        colorbar; view(2); shading interp; 
        caxis([ax_min, ax_max])
        title(sprintf('P%d_{2nd}', k))
        fprintf('P%d max / min: \t\t %f \t %f\n',  k, max(PN_2nd{(k+1)/2}.radiativeEnergy), min(PN_2nd{(k+1)/2}.radiativeEnergy))

        figure
        colormap(cmap);
        coarse2Fine = interpMeshFunction(PN_2nd{(k+1)/2}.radiativeEnergy, mesh, meshFine.points);
        trisurf(meshFine.connectivityList', meshFine.points(1, :), meshFine.points(2, :), coarse2Fine - PN_2ndFine{(k+1)/2}.radiativeEnergy );
        colorbar; view(2); shading interp
        title(sprintf('P%d_{2nd} - P%d_{2nd, fine}', k, k))
        fprintf('diffP%d max / min: \t %f \t %f\n', k, max(coarse2Fine - PN_2ndFine{(k+1)/2}.radiativeEnergy), min(coarse2Fine - PN_2ndFine{(k+1)/2}.radiativeEnergy))
    end

    %line plot
    figure
    st = {'b', 'r', 'g', 'k'};
    stIso = {'b--', 'r--', 'g--', 'y--'};
    for k = [1, 3, 5]
        plot(centerLine(1, :) , PN_2ndAlongLine{(k+1)/2}, st{(k+1)/2});
        hold on
        plot(centerLine(1, :) , PN_2ndIsoAlongLine{(k+1)/2}, stIso{(k+1)/2});
    end
    plotPwConstantFunction1D((DOAlongLine(1 : end - 1) + DOAlongLine(2 : end)) / 2, centerLine(1, :), 'k')
    
%     surface2tikz_plain(figure(1),'./testCase5_DO23','colorbar','horizontal','standalone',false,'dpi',500);
%     surface2tikz_plain(figure(3),'./testCase5_P1_2ndCoarse','colorbar','horizontal','standalone',false,'dpi',500);
%     surface2tikz_plain(figure(5),'./testCase5_P3_2ndCoarse','colorbar','horizontal','standalone',false,'dpi',500);
%     surface2tikz_plain(figure(7),'./testCase5_P5_2ndCoarse','colorbar','horizontal','standalone',false,'dpi',500);
end


function writeData2File()
    %% line plot
    fileName = [buildFolderMain, 'PN2ndAlongLine.csv'];
    headers = {'x', 'N1', 'N3', 'N5', 'N7', 'DO'};
    matrix = [centerLine(1, :)', PN_2ndAlongLine{1}', PN_2ndAlongLine{2}', PN_2ndAlongLine{3}', PN_2ndAlongLine{4}', DOAlongLine];
    write2File(fileName, headers, matrix)
    
    fileName = [buildFolderMain, 'PN2ndIsoAlongLine.csv'];
    headers = {'x', 'N1', 'N3', 'N5', 'N7'};
    matrix = [centerLine(1, :)', PN_2ndIsoAlongLine{1}', PN_2ndIsoAlongLine{2}', PN_2ndIsoAlongLine{3}', PN_2ndIsoAlongLine{4}'];
    write2File(fileName, headers, matrix)
    
    %% errors
    fileName = [buildFolder, 'errors.csv'];
    headers = {'fun', 'ToRefMesh', 'ToDO'};
    matrix = [[1, 3, 5, 7]', L2DifferenceRelFineMesh, L2DifferenceRelDO];
    fileID = fopen(fileName, 'w');
    fprintf(fileID, '%s,%s,%s\n', headers{1}, headers{2}, headers{3});
    for k = 1 : size(matrix, 1)
        tmp = sprintf('%d', matrix(k, 1));
        for j = 2 : size(matrix,  2)
           tmp = [tmp,sprintf('&\t') sprintf('%1.2e', matrix(k, j))];  
        end
        fprintf(fileID, [tmp, '\\\\\n']);
    end
    fclose(fileID);
    
end

end 

function [dPN2ndDOL2Rel] = relativeDistancePN2ndDO(PN_2nd, solutionDO)
    % eval coarse solution on fine mesh
    mesh.points = PN_2nd.points';
    mesh.connectivityList = double(PN_2nd.connectivityList' + 1);

    DONormL2 = L2NormElementFunction(solutionDO.radiativeEnergyDO, mesh.points, mesh.connectivityList);
    dPN2ndDOL2Rel = L2DifferenceMeshFunctionVsPwConstant(PN_2nd.radiativeEnergy, solutionDO.radiativeEnergyDO,...
                    mesh.points, mesh.connectivityList) / DONormL2;
end

function [dPN2ndRefInfRel, dPN2ndRefL2Rel] = relativeDistancePN2ndRefFine(fine, coarse)
    % eval coarse solution on fine mesh
    mesh.points = coarse.points';
    mesh.connectivityList = double(coarse.connectivityList' + 1);
    coarse2fine = interpMeshFunction(coarse.radiativeEnergy, mesh, fine.points');
    
    refNormInf = max(abs(fine.radiativeEnergy));
    dPNRefInfAbs = max(abs(coarse2fine - fine.radiativeEnergy));
    dPN2ndRefInfRel = dPNRefInfAbs  / refNormInf;
    
    meshFine.points = fine.points';
    meshFine.connectivityList = double(fine.connectivityList' + 1);
    refNormL2 = L2NormMeshFunction(fine.radiativeEnergy, meshFine.points, meshFine.connectivityList);
    dPN2ndRefL2Rel = L2NormMeshFunction(coarse2fine - fine.radiativeEnergy, meshFine.points, meshFine.connectivityList) ...
                      / refNormL2;
end

function buildFolder = pathDeclarations(nameTestCase)
    restoredefaultpath;
    buildFolder = ['./build/', nameTestCase, '/'];
    addpath('./testCases/')
    addpath(genpath('./SPN/')) 
    addpath(genpath('./discreteOrdinates/'))
    addpath('./tools/')
    addpath('./PN2nd/')
end


