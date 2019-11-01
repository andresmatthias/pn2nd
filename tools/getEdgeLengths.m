function edgeLengths = getEdgeLengths(points, edgeList)
% GETEDGELENGTHS Compute the lengths of edges in a triangular mesh.
% 
%   For a demonstration on how to use this function,
%   see also MESHAUXTEST
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%
P1 = points(:, edgeList(1, :));
P2 = points(:, edgeList(2, :));
edgeLengths = sqrt(sum((P2 - P1).^2, 1));
end
 