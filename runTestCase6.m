function runTestCase6()
% RUNTESTCASE6 Run the test case 6 and visualize the results.
% 
%   User instruction:
%       
%       1) Execute generateTestCase6GeoFile.m in folder
%           ./testCases/testCase6 to create the .geo file for this
%           customized domain.
%       2) Execute genMeshTestCase6.sh in folder ./testCases/testCase6 in
%           the command line (knowing dolfin and gmsh).
%       3) Use generateTestCase6.m first to generate second-order 
%           formulation of PN equations.
%       4) Run runTestCase5.py in Python (with FEniCS) in folder ./FEniCS .
%       5) Run runDiscreteOrdinatesTestCase5.m
%       6) Then execute this file.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

clear
close all
%% setup
nameTestCase = 'testCase6';
buildFolder = pathDeclarations(nameTestCase);
meshName = 'meshTestCase6_0';
meshNameFine = 'meshTestCase6_1';
maxExactDegree = [15, 23];

%% discrete ordinates
solutionDOCoarse = load(sprintf('%s%s_radEnergyDiscOrdinates_Deg_%d.mat', buildFolder, meshName, maxExactDegree(1)));
solutionDO = load(sprintf('%s%s_radEnergyDiscOrdinates_Deg_%d.mat', buildFolder, meshName, maxExactDegree(2)));

%% PN2nd
ax_min = min(solutionDO.radiativeEnergyDO);
ax_max = max(solutionDO.radiativeEnergyDO);
PN_2nd = {};
PN_2ndFine = {};
L2DifferenceRelFineMesh = zeros(4, 1);
L2DifferenceRelDO = zeros(4, 1);
L2DifferenceRelDODOCoarse = L2NormElementFunction(solutionDOCoarse.radiativeEnergyDO - solutionDO.radiativeEnergyDO, solutionDO.mesh.points, solutionDO.mesh.connectivityList) ...
    / L2NormElementFunction(solutionDO.radiativeEnergyDO, solutionDO.mesh.points, solutionDO.mesh.connectivityList);
fprintf('relative L2 difference DO: %1.2e\n', L2DifferenceRelDODOCoarse)

for k = [1, 3, 5, 7]
    PN_2nd{end + 1} = load(sprintf('%sFEniCS/radiativeEnergy_P%d2nd_%s.mat', buildFolder, k, meshName));
    PN_2ndFine{end + 1} = load(sprintf('%sFEniCS/radiativeEnergy_P%d2nd_%s.mat', buildFolder, k, meshNameFine));
    [dPN2ndRefInfRel, L2DifferenceRelFineMesh((k+1)/2)] = relativeDistancePN2ndRefFine(PN_2ndFine{end}, PN_2nd{end});
    L2DifferenceRelDO((k+1)/2) = relativeDistancePN2ndDO(PN_2nd{end}, solutionDO);
    fprintf('Relative difference (inf) P%d: ||coarse - fine||_inf / ||fine||_inf: \t%f\n', k, dPN2ndRefInfRel);
    fprintf('Relative difference (L2) P%d: ||coarse - fine||_L2 / ||fine||_L2: \t%f\n', k, L2DifferenceRelFineMesh((k+1)/2));
    fprintf('Relative difference (L2) ||P%d - DO||_L2 / ||DO||_L2: \t\t\t%f\n\n', k, L2DifferenceRelDO((k+1)/2));

    ax_min = min(ax_min, min(PN_2nd{end}.radiativeEnergy));
    ax_max = max(ax_max, max(PN_2nd{end}.radiativeEnergy));
end

          
%%
visualize();
writeErrors2File();

function writeErrors2File()
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

function visualize()
    cmap = viridis(256);

    mesh.points = PN_2nd{1}.points';
    mesh.connectivityList = double(PN_2nd{1}.connectivityList' + 1);
    meshFine.points = PN_2ndFine{1}.points';
    meshFine.connectivityList = double(PN_2ndFine{1}.connectivityList' + 1);
    fprintf('\nCoarse mesh: \t no points: %d, \t no elements: %d', size(mesh.points, 2), size(mesh.connectivityList, 2))
    fprintf('\nFine mesh: \t no points: %d, \t no elements: %d\n', size(meshFine.points, 2), size(meshFine.connectivityList, 2))

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
    
    
%     surface2tikz_plain(figure(1),'./testCase6_DO23','colorbar','horizontal','standalone',false,'dpi',500);
%     surface2tikz_plain(figure(3),'./testCase6_P1_2ndCoarse','colorbar','horizontal','standalone',false,'dpi',500);
%     surface2tikz_plain(figure(5),'./testCase6_P3_2ndCoarse','colorbar','horizontal','standalone',false,'dpi',500);
%     surface2tikz_plain(figure(7),'./testCase6_P5_2ndCoarse','colorbar','horizontal','standalone',false,'dpi',500);
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
    pathToTestCase = ['./testCases/', nameTestCase, '/'];
    addpath(pathToTestCase)
    addpath('./testCases/')
    addpath('./discreteOrdinates/2D/')
    addpath('./discreteOrdinates/')
    addpath('./tools/')
    addpath('./PN2nd/')
end