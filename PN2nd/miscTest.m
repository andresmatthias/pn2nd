% MISCTEST Perform miscellaneous unit tests in the context of the 
%   second-order formulation of the PN equations.
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

%% own cholesky solver
n = 50;
A = randn(n, n);
b = randn(n, 1);
b = A' * b;
A = A' * A;
xMatlab = A \ b;
L = chol(A, 'lower');
yMatlabChol = L \ b;
xMatlabChol = L' \ yMatlabChol;
x = choleskySolve(A, b);
diff1 = norm(x - xMatlabChol, 'inf');
diff2 = norm(x - xMatlab, 'inf');
if (diff1 >= 1e-10) || (diff2 >= 1e-10)
    fprintf('Difference custom cholesky vs. Matlab cholesky: %f\nDifference custom cholesky vs. Matlab backslash: %f', diff1, diff2)
end
assert(diff1 < 1e-10)
assert(diff2 < 1e-10)

%% rotation matrix
addpath('quadratureRules/')
for k = 1 : 50
    a = randn(3, 1);
    a = a / norm(a);
    b = randn(3, 1);
    b = b / norm(b);
    R = getGeomRotationa2b(a, b);
    assert(norm(R * a - b) < 1e-12)
    assert(abs(det(R) - 1) < 1e-12)
end

%% parametrization of angular variable Omega
n = 50;
Omega = randn(3, n);
Omega = Omega ./ sqrt(sum(Omega.^2, 1));
[mu, phi] = sphereCartesian2Param(Omega);
Omega1 = sphereParam2Cartesian(mu, phi);
assert(norm(Omega1 - Omega, 'inf') < 1e-14)

mu = (rand(1, n) - 1 / 2) * 2;
phi = rand(1, n) * 2 * pi;
Omega = sphereParam2Cartesian(mu, phi);
[mu1, phi1] = sphereCartesian2Param(Omega);
assert(max(abs(mu - mu1)) < 1e-14)
phi1(phi1 < 0) = phi1(phi1 < 0) + 2 * pi; %[-pi, pi] -> [0, 2 * pi]
phi1(abs(phi1 - 2 * pi) < 1e-15) = 0;
phi(abs(phi - 2 * pi) < 1e-15) = 0;
assert(max(abs(phi - phi1)) < 1e-14)

%% reflection at surface
outerNormalVectors3D = [1, 0, 0;
                        0, 1, 0;
                        -1, 0, 0;
                        0, -1, 0;
                        0, 0, 1;
                        0, 0, -1]';
Omega = [-1/sqrt(2), 1/sqrt(2), 0;
         -1/sqrt(2), -1/sqrt(2), 0;
          1/sqrt(2), -1/sqrt(2), 0;
          1/sqrt(2), 1/sqrt(2), 0;
          1/sqrt(2), 0, -1/sqrt(2);
          1/sqrt(2), 0, 1/sqrt(2)]';                   
                   
OmegaRefTrue = [1/sqrt(2), 1/sqrt(2), 0;
                -1/sqrt(2), 1/sqrt(2), 0;
                -1/sqrt(2), -1/sqrt(2), 0;
                1/sqrt(2), -1/sqrt(2), 0;
                1/sqrt(2), 0, 1/sqrt(2);
                1/sqrt(2), 0, -1/sqrt(2);]';

OmegaRef = reflectionAtSurface(outerNormalVectors3D, Omega);
assert(norm(OmegaRef - OmegaRefTrue, 'inf') < 1e-15)
            

%% integral outer product
x = sym('x');
f = [x; x.^2];
g = [x.^3; x.^4];
outerP = f * g.';
ITrue = double(int(outerP, x, -1, 1));
GaussQuadX = [-sqrt(3 / 7 + 2 / 7 * sqrt(6 / 5));
              -sqrt(3 / 7 - 2 / 7 * sqrt(6 / 5));
               sqrt(3 / 7 - 2 / 7 * sqrt(6 / 5));
               sqrt(3 / 7 + 2 / 7 * sqrt(6 / 5))]';
GaussQuadWeights = [(18 - sqrt(30)) / 36;
                    (18 + sqrt(30)) / 36;
                    (18 + sqrt(30)) / 36;
                    (18 - sqrt(30)) / 36;]';
f = matlabFunction(f);
g = matlabFunction(g);
fAtQuad = f(GaussQuadX);
gAtQuad = g(GaussQuadX);
I = integralOuterProduct(fAtQuad, gAtQuad, GaussQuadWeights);
assert(norm(I - ITrue, 'inf') < 1e-14);
                    
                    
%% projection of scattering operator
addpath('./basisFunctions/')
for spatialDimension = [1, 2, 3]
    N = 10;
    % this is not a scattering kernel, but we use it to validate computation of
    % the double integral
    kernelForm = @(vx, vy, vz, wx, wy, wz) ones(size(vx .* wx)); 
    CScatter = projectScatteringOperator(N, kernelForm, spatialDimension);
    CScatterTrue = zeros(size(CScatter));
    CScatterTrue(1,1) = 4 * pi;
    assert(norm(CScatter - CScatterTrue, 'inf') < 1e-12)
end             
                    
                    