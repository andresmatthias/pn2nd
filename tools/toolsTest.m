% TOOLSTEST Perform miscellaneous unit tests for different tools 
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

%% L2 norm of piecewise-constant function on triangular mesh
points = [0, 0; 1, 0; 1, 1; 0, 1]';
connectivityList = [1, 2, 4; 2, 4, 3]';
f = [ -2 / 3, 4 / 3];
E = L2NormElementFunction(f, points, connectivityList);
intTrue = sqrt( 1/2 * (-2/3)^2 + 1/2 * (4/3)^2);
assert(abs(E - intTrue) < 1e-15)

%% L2 norm piecewise-linear function on triangular mesh
points = [0, 0; 1, 0; 1, 1; 0, 1]';
connectivityList = [1, 2, 4; 2, 4, 3]';
syms x y
f1 = -1 + 2*x + y;
f2 = -1 + 2*x + y;
intTrue = double( int(int(f1^2, y, 0, 1 - x), x, 0, 1) ...
                 +int(int(f2^2, y, 1 - x, 1), x, 0, 1));

f = [-1, 1, 2, 0];
E = L2NormMeshFunction(f, points, connectivityList);
assert(abs(E - sqrt(intTrue)) < 1e-15)

%% L2 difference piecewise-linear vs. piecewise const.
points = [0, 0; 1, 0; 1, 1; 0, 1]';
connectivityList = [1, 2, 4; 2, 4, 3]';
syms x y
f1 = -1 + 2*x + y;
f2 = -1 + 2*x + y;

g1 = 3;
g2 = 4;

diffTrue = double(int(int((f1 - g1)^2, y, 0, 1 - x), x, 0, 1) ...
                 +int(int((f2 - g2)^2, y, 1 - x, 1), x, 0, 1));

f = [-1, 1, 2, 0];
g = [3, 4];
E = L2DifferenceMeshFunctionVsPwConstant(f, g, points, connectivityList);
assert(abs(E - sqrt(diffTrue)) < 1e-15)