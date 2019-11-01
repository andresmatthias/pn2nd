% DISCRETEORDINATES2DTEST Perform unit tests regarding the discrete
%   ordinates method in 2D.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%
%% geodesic edge lengths
A = rand(3, 100);
A = A ./ sqrt(sum(A.^2, 1)); 
B = rand(3, 100);
B = B ./ sqrt(sum(B.^2, 1)); 

arcLengthExpected = acos(sum(A .* B, 1));
arcLength = geodesicEdgeLengthsUnitSphere(A, B);
assert(norm(arcLength - arcLengthExpected, 'inf') < 1e-13);

%% distance based interpolation unit sphere
pointsOnSphereCartesian = [1, 0, 0;
                           0, 1, 0;
                           0, 0, 1]';
queryPointsCartesian = [ [1, 1, 1]'/sqrt(3), [1, 0, 0]'];
[~, weights] = distanceBasedInterpolationUnitSphere(pointsOnSphereCartesian, queryPointsCartesian);
assert(norm(weights(:, 1) - [1, 1, 1]' / 3) < 1e-15)
assert(norm(weights(:, 2) - [1, 0, 0]') < 1e-15)

%% weights interpolation
pointsOnSphereCartesian = [1, 0, 0;
                           0, 1, 0;
                           0, 0, 1]';
queryPointsCartesian = rand(3, 100);
queryPointsCartesian = queryPointsCartesian ./ sqrt(sum(queryPointsCartesian.^2, 1)); 
weights = distanceBasedInterpolationUnitSphere(pointsOnSphereCartesian, queryPointsCartesian);

%% area of triangles on unit sphere
A = rand(3, 100);
A = A ./ sqrt(sum(A.^2, 1)); 
B = rand(3, 100);
B = B ./ sqrt(sum(B.^2, 1)); 
C = rand(3, 100);
C = C ./ sqrt(sum(C.^2, 1)); 

area = areaUnitSphereTriangle(A, B, C);
solidAngle = solidAngleUnitSphereTriangle(A, B, C);
assert(norm(solidAngle - area, 'inf') < 1e-12)

%% solid angle for validation
A = [1, 0, 0]';
B = [0, 1, 0]'; 
C = [0, 0, 1]';
solidAngle = solidAngleUnitSphereTriangle(A, B, C);
assert(abs(solidAngle - 4 * pi / 8) < 1e-14)

