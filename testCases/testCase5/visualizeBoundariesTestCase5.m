function visualizeBoundariesTestCase5()
% VISUALIZEBOUNDARIESTESTCASE5 Visualize the boundary IDs for test case 5.
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

mesh = load('../../build/testCase5/meshTestCase5_0.mat');
aux(mesh.connectivityList, mesh.points, mesh.edgeList, mesh.boundaryEdgesIdx2Id)

end

function [] = aux(connectivityList, points, edgeList,...
    boundaryEdgesIdx2ID)
fprintf('\nVisualize boundaries\n')
close all
trimesh(connectivityList', points(1, :)', points(2, :)','Color','k')
hold on
sp = zeros(4, 1);
for k = 1:length(boundaryEdgesIdx2ID.edgeIdx)
    p = points(:, edgeList(:, boundaryEdgesIdx2ID.edgeIdx(k)));
    switch boundaryEdgesIdx2ID.boundaryId(k)
        case 1
            sp(1) = plot(p(1, :), p(2, :), 'r', 'LineWidth', 3);
        case 2
            sp(2) = plot(p(1, :), p(2, :), 'b', 'LineWidth', 3);
        case 3
            sp(3) = plot(p(1, :), p(2, :), 'g', 'LineWidth', 3);
        case 4
            sp(4) = plot(p(1, :), p(2, :), 'y', 'LineWidth', 3);
    end
end
legend(sp, '1','2','3','4')
end
