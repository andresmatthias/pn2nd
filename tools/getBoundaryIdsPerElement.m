function boundaryIdsPerElement = getBoundaryIdsPerElement(nr_elements, ...
    edgesPerElement, boundaryEdgesIdx2Id, elementNeighbors)
% GETBOUNDARYIDSPERELEMENT Assign boundary IDs to the edges of an element.
% 
%   For each element, the k-th (k=1,2,3) edge is opposite of the k-th point
%   in the connectivity list. Assign the ID of the corresponding boundary,
%   if the edge is at the boundary, and NaN otherwise. 
% 
%   E.g., T_k = [A, B, C], edge AB is boundary edge with boundaryID XXX,
%   the other edges are interior edges, then we get [NaN;NaN;XXX].
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

fprintf('\nGet boundary edges per element ...\n')
boundaryIdsPerElement = nan(3, nr_elements);
for k = 1:nr_elements
    edges_k = edgesPerElement(:, k);
    [idxa, idxb] = ismember(edges_k, boundaryEdgesIdx2Id.edgeIdx);
    idxb = idxb(idxa); % to kill zero indices
    boundaryIdsPerElement(idxa, k) = boundaryEdgesIdx2Id.boundaryId(idxb);
end

checksum = length(find(~isnan(boundaryIdsPerElement(1, :))))...
    +length(find(~isnan(boundaryIdsPerElement(2, :))))...
    +length(find(~isnan(boundaryIdsPerElement(3, :))));
if checksum ~= length(boundaryEdgesIdx2Id.edgeIdx)
    error('Something went wrong with matching the boundary edges!')
end

tmp = elementNeighbors + boundaryIdsPerElement;
if any(~isnan(tmp))
    % if boundary edge: no adjacent element (-> nan)
    % if interior element: no boundary edge (-> nan)
    error('Boundary edges not aligned with boundary elements!')
end
end