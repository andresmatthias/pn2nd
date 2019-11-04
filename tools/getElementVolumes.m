function elementVol = getElementVolumes(points, connectivityList)
% GETELEMENTVOLUMES Get the volumes (2D) of the elements of a triangular
%   mesh.
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
A = points(:, connectivityList(1, :));
B = points(:, connectivityList(2, :));
C = points(:, connectivityList(3, :));
eAB = B-A;
eAC = C-A;
elementVol = abs( eAB(1, :) .* eAC(2, :) - eAB(2, :) .* eAC(1, :)) / 2;
end