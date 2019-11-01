function CScatter = projectScatteringOperator(N, kernelFun, spatialDimension, varargin)
% PROJECTSCATTERINGOPERATOR Compute the moments of the scattering operator
%   w.r.t. basis functions (real spherical harmonics) up to degree N. 
%   Choose a subset of basis functions according to symmetry assumptions 
%   (spatial dimension).
%
%   The kernel function needs to allow vectorized evaluation, e.g., 
%   k = @(vx, vy, vz, wx, wy, wz) 1 / 4 / pi * ones(size(vx .* wx))
%
% See also MISCTEST, PN2NDREDUCTIONTEST
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

%% choose quadrature for scattering operator
% We need to integrate function of degree deg(b)*2 + deg(kernel) (deg(b) = N)
% addpath('./quadratureRules/')
% addpath('./basisFunctions/')
if nargin > 3
    degKernel = varargin{1};
else
    degKernel = 50; % not known in advance, so take this value for now
end
fprintf('\nassume degree(scatter_kernel) <= %d\n', degKernel);
maxExactDegree = degKernel + 2 * N;
[weights, mu, phi] = sphericalQuadratureFull(maxExactDegree);

basisAtQuad = evalONB([mu; phi], N, spatialDimension);
nMoments = getNumberOfBasisFunctions(N, spatialDimension);
Omega = sphereParam2Cartesian(mu, phi);
wx = Omega(1, :); % row vector
wy = Omega(2, :);
wz = Omega(3, :);
vx = wx'; vy = wy'; vz = wz'; % column vector
kEval = kernelFun(vx, vy, vz, wx, wy, wz); % matrix
basisTimesWeights = repmat(weights, nMoments, 1) .* basisAtQuad;
CScatter = basisTimesWeights * (kEval * basisTimesWeights');

%% truncate
% The weights of the quadrature rule sum up to 4*pi + small_value ; from 
% this we cannot assume exact computation of zeros. Rounding the computed 
% integrals improves the sparsity-structure of the matrix.

exponent10 = log(abs(sum(weights) - 4 * pi)) / log(10);
cut = min(16, floor(-exponent10) - 1);
fprintf('Cut off scattering integrals after %d digits\n', cut)
CScatter = round(CScatter, cut);  % let 'almost zeros' caused by quadrature disappear
end