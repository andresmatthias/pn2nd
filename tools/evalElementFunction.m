function [fAtQueryPoints] = evalElementFunction(f, mesh, queryPoints)
% EVALELEMENTFUNCTION Evaluate a function, piecewise constant w.r.t. 
%   elements of a triangular mesh, at given query points.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
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