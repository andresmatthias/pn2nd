function visualizeBoundariesTestCase6()
% VISUALIZEBOUNDARIESTESTCASE6 Visualize the boundary IDs for test case 6.
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

mesh = load('../../build/testCase6/meshTestCase6_0.mat');
aux(mesh.connectivityList, mesh.points, mesh.edgeList, mesh.boundaryEdgesIdx2Id)
end

function [] = aux(connectivityList, points, edgeList,...
    boundaryEdgesIdx2ID)
fprintf('\nVisualize boundaries\n')
close all
trimesh(connectivityList', points(1, :)', points(2, :)','Color','k')
hold on
sp = cell(4, 1);

cmap = parula(length(unique(boundaryEdgesIdx2ID.boundaryId)));

for k = 1:length(boundaryEdgesIdx2ID.edgeIdx)
    
    p = points(:, edgeList(:, boundaryEdgesIdx2ID.edgeIdx(k)));
    idx = boundaryEdgesIdx2ID.boundaryId(k);
    c = cmap(idx, :);
    sp{idx} = plot(p(1, :), p(2, :), 'color', c, 'LineWidth', 3);
end

legendAux1 = [sp{1}];
legendAux2 = {sprintf('%d', 1)};
for k = 2 : size(cmap, 1)
   legendAux1 = [legendAux1, sp{k}]; 
   legendAux2{end+1} = sprintf('%d', k);
end
legend(legendAux1, legendAux2)
end
