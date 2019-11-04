function [fAtQueryPoints] = evalElementFunction(f, mesh, queryPoints)
% EVALELEMENTFUNCTION Evaluate a function, piecewise constant w.r.t. 
%   elements of a triangular mesh, at given query points.
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%
if size(mesh.connectivityList, 1) == 3 && size(mesh.points, 1) == 2 ...
    && size(queryPoints, 1) == 2 && isvector(f) 
    TR = triangulation(mesh.connectivityList', mesh.points');
    [idx, ~] = pointLocation(TR, queryPoints');
    idxNan = find(isnan(idx));
    for k = 1 : length(idxNan)
        while isnan(idx(idxNan(k)))
            fprintf('disturb point\n')
            [idxAux, ~] = pointLocation(TR, queryPoints(:, idxNan(k))' + 1e-14 * randn(1, 2));
            idx(idxNan(k)) = idxAux;
        end
    end

    fAtQueryPoints = f(idx);
else
    error('Dimension error in connectivity list!') 
end
end