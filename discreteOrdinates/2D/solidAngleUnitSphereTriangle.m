function solidAngle = solidAngleUnitSphereTriangle(A, B, C)
% SOLIDANGLEUNITSPHERETRIANGLE Compute the solid angle of a geodesic
%   triangle on the unit sphere (only for small angles).
% 
%   Only for small angles, e.g. we get an error for 
%   A = [1, 0 , 0]', B = [-1, 0, 0]', C = [0, 0, 1]'
%
%   For a demonstration on how to use this function,
%   see also DISCRETEORDINATES2DTEST
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%
s = abs(dot(A, cross(B, C)));
a = sqrt(sum(A.^2, 1));
b = sqrt(sum(B.^2, 1));
c = sqrt(sum(C.^2, 1));
solidAngle = 2 * atan(s ./ (a .* b .* c + dot(A, B) .* c + dot(A, C) .* b + dot(B, C) .* a));
end