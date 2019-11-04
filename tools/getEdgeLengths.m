function edgeLengths = getEdgeLengths(points, edgeList)
% GETEDGELENGTHS Compute the lengths of edges in a triangular mesh.
% 
%   For a demonstration on how to use this function,
%   see also MESHAUXTEST
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%
P1 = points(:, edgeList(1, :));
P2 = points(:, edgeList(2, :));
edgeLengths = sqrt(sum((P2 - P1).^2, 1));
end
 