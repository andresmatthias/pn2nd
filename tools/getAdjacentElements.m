function elementNeighbors = getAdjacentElements(triang)
% GETADJACENTELEMENTS For a given triangulation of a 2D mesh, get the
%   indices of the adjacent elements. 
% 
%   For an element (cell), the k-th (k = 1,2,3) adjacent element is located
%   opposite the k-th point in the connectivity list (analog to the order
%   of the edges per element.)
%   NaN is assigned in case there is no adjacent element (at the boundary 
%   of the domain).
% 
%   For a demonstration on how to use this function,
%   see also MESHAUXTEST
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%
%% adjacent elements
fprintf('\nGet adjacent elements ...\n')
elementNeighbors = neighbors(triang);
elementNeighbors = elementNeighbors';
end