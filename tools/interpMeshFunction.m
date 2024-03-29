function [fAtQueryPoints] = interpMeshFunction(f, mesh, queryPoints)
% INTERPMESHFUNCTION Evaluate a piecewise linear function w.r.t. a
%   triangular mesh, given by its nodal values, at query points by
%   barycentric interpolation (without extrapolation, i.e., query points
%   have to be inside the domain).
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
    [idx, weights] = pointLocation(TR, queryPoints');
    idxNan = find(isnan(idx));
    for k = 1 : length(idxNan)
        while isnan(idx(idxNan(k)))
            [idxAux, weightsAux] = pointLocation(TR, queryPoints(:, idxNan(k))' + 1e-14 * randn(1, 2));
            idx(idxNan(k)) = idxAux;
            weights(idxNan(k), :) = weightsAux;
        end
    end

    fAtQueryPoints = sum(weights' .* f(mesh.connectivityList(:, idx)), 1);
else
    error('Dimension error in connectivity list!') 
end
end