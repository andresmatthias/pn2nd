function generateGmshGeoFile(boundaryPolygon, meshName, lc, geoGenParams)
% GENERATEGMSHGEOFILE Generate a .geo file for Gmsh, based on a boundary
%   polygon and characteristic length.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%
fileID = fopen(sprintf('./%s.geo', meshName), 'w');
fprintf(fileID, '// 2D geometry of testCase2D_4\n\n');

if norm(boundaryPolygon(1, :) - boundaryPolygon(end, :)) > 1e-14
   error('First and last point must coincide!') 
end

nBoundary = size(boundaryPolygon, 1) - 1;
%% points
for k = 1 : size(boundaryPolygon, 1)
   fprintf(fileID, 'Point(%d) = {%f, %f, 0, %f};\n', k, boundaryPolygon(k, 1), boundaryPolygon(k, 2), lc(k)); 
end
fprintf(fileID, '\n');
%% lines
for k = 1 : nBoundary
   if k < nBoundary
       next = k + 1;
   else
       next = 1;
   end
   fprintf(fileID, 'Line(%d) = {%d, %d};\n', k, k, next); 
end
fprintf(fileID, '\n');
%%
aux = '1';
for k = 2 : nBoundary
   aux = [aux, ', ', int2str(k)]; 
end
fprintf(fileID, 'Line Loop(1) = {%s};\n\n', aux);
fprintf(fileID, 'Plane Surface(1) = {1};\n\n');
fprintf(fileID, 'Physical Point(1) = {%s};\n', aux);

for k = 1 : nBoundary
   fprintf(fileID, 'Physical Line(%d) = {%d};\n', k, k);
end

fprintf(fileID, '\nPhysical Surface("2DLayer", 999) = {1};\n');

fclose(fileID);

%%
outerNormalVectors3D = zeros(3, nBoundary);
for k = 1 : nBoundary
    if k < nBoundary
       next = k + 1;
   else
       next = 1;
   end
    edge = boundaryPolygon(next,:) - boundaryPolygon(k,:);
    edge = edge / norm(edge);
    outerNormalVectors3D(:, k) = [rightNormal(edge');0];
end
boundaryIDs = (1 : nBoundary)';
% system(sprintf('gmsh %s.geo -2 -o %s.msh', meshname, meshname));
save(sprintf('%sNormals', meshName), 'outerNormalVectors3D', 'boundaryIDs', 'geoGenParams');
end

function p = rightNormal(pIn)
    p = [pIn(2, :);
         -pIn(1, :)];
end