% BASISFUNCTIONSTEST Perform unit tests regarding properties of the real 
%   spherical harmonics as orthonormal basis on the unit sphere, depending
%   on symmetry assumptions (spatialDimension)).
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

%% basis is orthonormal
addpath('../quadratureRules/')
addpath('../')
N = 10;
maxExactDegree = 2 * N; % deg(basis_i) <= N; consider basis_i * basis_j
for spatialDimension = [1, 2, 3]
    [weights, mu, phi] = sphericalQuadratureFull(maxExactDegree);
    b = evalONB([mu; phi], N, spatialDimension);
    nMoments = getNumberOfBasisFunctions(N, spatialDimension);
    I = integralOuterProduct(b, b, weights); % is of degree 2*N
    assert(max(max( abs(I -  eye(nMoments)) )) < 1e-13)
end 

%% even and odd linear indices are partition of all linear indices
N = 99;
for spatialDimension = [1, 2, 3]
   nMoments = getNumberOfBasisFunctions(N, spatialDimension);
   linIdx = 0 : nMoments - 1;
   idxEven = linearIdxOfEvenBasis(N, spatialDimension);
   idxOdd = linearIdxOfOddBasis(N, spatialDimension);
   assert(norm(idxEven - setdiff(linIdx, idxOdd)) == 0)
end

%% associated Legendre polynomials correct
% factorial in symbolic expression leads to strange results for large N (N=7)
N = 5;
mu = sym('mu');
phi = sym('phi');
P = evalLegendreAssociatedPolynomials([mu; phi], N);

diffMax = 0;
Pexact = cell(N + 1);

for l = 0 : N    
    for m = 0 : l
        % Without Condon-Shortley phase!!!
        Pexact{l + 1, m + 1} = simplify(1 / 2^l / factorial(l) * (1 - mu^2)^(m / 2) * diff((mu^2 - 1)^l, mu, l + m));
        assert( isequal(simplify(Pexact{l + 1, m + 1} - P{l + 1, m + 1}), sym(0)))      
    end
end

%% restriction of 3D spherical harmonics yields legendre polynomials
% evaluate for phi = 0 and compare to legendre polynomial
N = 20;
spatialDimension = 3;
pointsOnSphere = [linspace(-1, 1, 100); zeros(1, 100)];
sphHarmAtPoints = evalONB(pointsOnSphere, N, spatialDimension);
linIdx = 0 : getNumberOfBasisFunctions(N, spatialDimension) - 1;
[l, m] = linearIdx2DegOrder(linIdx, spatialDimension);
P = evalLegendreAssociatedPolynomials(pointsOnSphere(1, :), N);
differenceMax = zeros(length(l), 1);
for i = 1:length(linIdx)
   phi = 0;
   Y_desired = @(l, m) sqrt((2 * l + 1) *  factorial(l - m) / (4 * pi * factorial(l + m))) * P{l + 1, m + 1};
   if m(i) > 0
       S_desired = Y_desired(l(i), m(i)) * sqrt(2) * cos(m(i) * phi);
   elseif m(i) == 0
       S_desired = Y_desired(l(i), m(i));
   elseif m(i) < 0
       S_desired = Y_desired(l(i), abs(m(i)))*sqrt(2) * sin(abs(m(i)) * phi);
   end
   differenceMax(i) = max(abs(S_desired - sphHarmAtPoints(i, :)));
end
assert(max(differenceMax) < 1e-15)

%% correct number of basis functions (after reduction)
N = 20;
for spatialDimension = [1, 2, 3]
    pointsOnSphere = [linspace(-1, 1, 100); zeros(1, 100)];
    basisAtPoints = evalONB(pointsOnSphere, N, spatialDimension);
    assert(isequal(size(basisAtPoints, 1), getNumberOfBasisFunctions(N, spatialDimension)))
end