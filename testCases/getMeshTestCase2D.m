function [mesh] = getMeshTestCase2D(folderGmshFile, meshName, buildFolder,...
    getBoundaryEdgesIdx2IdTestCase2D_x)
[points, connectivityList] = importGmsh([folderGmshFile, meshName, '.msh'], '2D');

%% elements
[connectivityList] = orderConnectivityListCcw(points, connectivityList); % needed for normals

TR = triangulation(connectivityList', points');

mesh.name = meshName;
mesh.nElements = size(connectivityList, 2);
mesh.points = points;
mesh.connectivityList = connectivityList;
mesh.elementVolume = getElementVolumes(points, connectivityList);
mesh.edgeList = (sort(edges(TR), 2))';
mesh.edgeLengths = getEdgeLengths(points, mesh.edgeList);
mesh.edgesPerElement = getEdgesPerElement(connectivityList, mesh.edgeList);
mesh.elementNeighbors = getAdjacentElements(TR);
mesh.edgeNormalsPerElement = getEdgeNormalsPerElement(connectivityList, points);
mesh.boundaryEdgesIdx2Id = getBoundaryEdgesIdx2IdTestCase2D_x(TR, points, mesh.edgeList);
mesh.boundaryIdsPerElement = getBoundaryIdsPerElement(mesh.nElements, ...
    mesh.edgesPerElement, mesh.boundaryEdgesIdx2Id, mesh.elementNeighbors);

if ~isempty(buildFolder)
    save([buildFolder, meshName, '.mat'], '-struct', 'mesh')
end
end


