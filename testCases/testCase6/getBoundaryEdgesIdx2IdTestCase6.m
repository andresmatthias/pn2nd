function boundaryEdges = getBoundaryEdgesIdx2IdTestCase6(TR, points, edgeList, geoGenParams)
% GETBOUNDARYEDGESIDX2IDTESTCASE6 Map the index of an edge within its edge
%   list to the id of the boundary it belongs to.  
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

tol = 1e-7;
boundaryEdgesByPoints = freeBoundary(TR);
boundaryEdgesByPoints = (sort(boundaryEdgesByPoints, 2))'; % to use ismember(,'row') later
[idxa, boundaryEdgeIdx] = ismember(boundaryEdgesByPoints', edgeList', 'rows');
if any(idxa==0)
    error('Boundary edge not found!')
end

nrEdges = size(boundaryEdgeIdx, 1);
boundaryEdgeId = zeros(nrEdges, 1);
% solve this more elegant
P1 = points(:, edgeList(1, boundaryEdgeIdx));
P2 = points(:, edgeList(2, boundaryEdgeIdx));
midpoint = (P1 + P2) / 2;

%% copy from geo-file generator

boundaryMidpoints = (geoGenParams.boundaryPolygon(1 : end-1, :) + geoGenParams.boundaryPolygon(2 : end, :)) / 2;
boundaryMidpoints = boundaryMidpoints';

for k  = 1 : nrEdges
   [~, idx] = min(sum((boundaryMidpoints - midpoint(:, k)).^2, 1)); 
   boundaryEdgeId(k) = idx;
end

boundaryEdges.edgeIdx = boundaryEdgeIdx';
boundaryEdges.boundaryId = boundaryEdgeId;
end