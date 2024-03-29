function edgeNormalsPerElement = getEdgeNormalsPerElement(connectivityList, points)
% GETEDGENORMALSPERELEMENT Compute the outer normals per element of a
%   triangular mesh.
% 
%   For each element, the k-th (k=1,2,3) outer normal vector corresponds to
%   the edge opposite of the k-th point in the connectivity list.
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

fprintf('\nGet normals ...\n')
nElements = size(connectivityList, 2);
edgeNormalsPerElement = zeros(2, 3, nElements);
for k = 1 : nElements
    points_k = points(:, connectivityList(:, k));
    % this works only if connectivity_list is ordered ccw!!
    nAB = getNormalRight(points_k(:, 2) - points_k(:, 1));
    nBC = getNormalRight(points_k(:, 3) - points_k(:, 2));
    nCA = getNormalRight(points_k(:, 1) - points_k(:, 3));
    % opposite of points A,B,C
    edgeNormalsPerElement(:, 1, k) = nBC;
    edgeNormalsPerElement(:, 2, k) = nCA;
    edgeNormalsPerElement(:, 3, k) = nAB;
end
end
