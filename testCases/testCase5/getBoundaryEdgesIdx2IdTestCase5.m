function boundaryEdges = getBoundaryEdgesIdx2IdTestCase5(TR, points, edgeList)
% GETBOUNDARYEDGESIDX2IDTESTCASE5 Map the index of an edge within its edge
%   list to the id of the boundary it belongs to.  
%
%   The Gmsh geometry (rectangle) has the following boundary IDs:
%       1: right
%       2: top
%       3: left
%       4: bottom
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

nrEdges = size(boundaryEdgeIdx, 2);
boundaryEdgeId = zeros(nrEdges, 1);
% solve this more elegant
P1 = points(:, edgeList(1, boundaryEdgeIdx));
P2 = points(:, edgeList(2, boundaryEdgeIdx));
midpoint = (P1 + P2) / 2;
xmin = min(P1(1, :));
xmax = max(P1(1, :));
ymin = min(P1(2, :));
ymax = max(P1(2, :));
boundaryEdgeId(abs(midpoint(1, :) - xmin) < tol) = 3;
boundaryEdgeId(abs(midpoint(1, :) - xmax) < tol) = 1;
boundaryEdgeId(abs(midpoint(2, :) - ymin) < tol) = 4;
boundaryEdgeId(abs(midpoint(2, :) - ymax) < tol) = 2;
if any(boundaryEdgeId == 0)
    error('No matching boundary ID!')
end

boundaryEdges.edgeIdx = boundaryEdgeIdx';
boundaryEdges.boundaryId = boundaryEdgeId;
end