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
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%
%% adjacent elements
fprintf('\nGet adjacent elements ...\n')
elementNeighbors = neighbors(triang);
elementNeighbors = elementNeighbors';
end