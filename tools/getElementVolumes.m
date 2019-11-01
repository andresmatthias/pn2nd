function elementVol = getElementVolumes(points, connectivityList)
% GETELEMENTVOLUMES Get the volumes (2D) of the elements of a triangular
%   mesh.
% 
%   For a demonstration on how to use this function,
%   see also MESHAUXTEST
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%
A = points(:, connectivityList(1, :));
B = points(:, connectivityList(2, :));
C = points(:, connectivityList(3, :));
eAB = B-A;
eAC = C-A;
elementVol = abs( eAB(1, :) .* eAC(2, :) - eAB(2, :) .* eAC(1, :)) / 2;
end