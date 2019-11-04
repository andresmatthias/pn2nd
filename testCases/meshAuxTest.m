% MESHAUXTEST Perform unit tests regarding auxiliary structures of 2D
%   meshs.
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

%% counter-clockwise ordering 2D mesh
addpath('../tools/')
points = [0, 1, 0;
          0, 0, 1];
connectivityList = [1, 2, 3]';
connectivityList_clockwise = [1, 3, 2]';
aux = orderConnectivityListCcw(points, connectivityList);
assert(all(connectivityList == aux))
aux = orderConnectivityListCcw(points, connectivityList_clockwise);
assert(all(connectivityList == aux))

%% enhance structure 2D mesh
addpath('../tools/')
% element volumes 2D mesh
[points, connectivityList] = importGmsh(['./', 'testMesh', '.msh'], '2D');
TR = triangulation(connectivityList', points');
elementVolumes = getElementVolumes(points, connectivityList);
assert(abs(sum(elementVolumes) - 1) < 1e-15);

% edge lengths in 2D mesh
edgeList = (sort(edges(TR), 2))';
edgeLengths =  getEdgeLengths(points, edgeList);
edgeLengthsSample = getEdgeLengths(points, [1, 2; 1, 5]');
assert(norm(edgeLengthsSample - [1, sqrt(2) / 2]) < 1e-15);

% edges per element in 2D mesh
% check via Heron, if edge lengths and element size correspond
% https://de.wikipedia.org/wiki/Dreiecksfl%C3%A4che
edgesPerElement = getEdgesPerElement(connectivityList, edgeList);
s = sum(edgeLengths(edgesPerElement), 1) / 2;
F = sqrt(s .* (s - edgeLengths(edgesPerElement(1, :))) ...
    .* (s - edgeLengths(edgesPerElement(2, :))) ...
    .* (s - edgeLengths(edgesPerElement(3, :))));
assert(norm(elementVolumes - F, 'inf') < 1e-13)

% adajent elements in 2D mesh
% check if order in element_edges is aligned with order in
% element_neighbors
elementNeighbors = getAdjacentElements(TR);
nElements = size(connectivityList, 2);
for  k = 1 : nElements
    for  j = 1:3
        edge_k_by_points = edgeList(:, edgesPerElement(j, k));
        neighborElementIdx = elementNeighbors(j, k);
        if neighborElementIdx > 0 % check for boundary element
            pointsNeighbor = connectivityList(:, neighborElementIdx);
            tmp = [sort(pointsNeighbor([1,2]))'; sort(pointsNeighbor([2,3]))'; sort(pointsNeighbor([1,3]))'];
            assert(ismember(edge_k_by_points', tmp, 'rows'))
        end
    end
end

% edge normals in 2D mesh
edgeNormalsPerElement = getEdgeNormalsPerElement(connectivityList, points);
for k = 1 : nElements
    for j = 1 : 3
       A = points(:, edgeList(1, edgesPerElement(j, k))); 
       B = points(:, edgeList(2, edgesPerElement(j, k)));
       assert(abs(dot(B - A, edgeNormalsPerElement(:, j, k))) < 1e-15);
    end
end