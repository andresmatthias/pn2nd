function [realSphericalHarmonicsAtPoints] = evalRealSphericalHarmonics(pointsOnSphereMuPhi, lmax)
% EVALREALSPHERICALHARMONICS Evaluate all real spherical harmonics up to 
%   degree lmax at given points on sphere.
%
%   In our implementation we omit the Condon-Shortley phase factor in
%   the definition of Legendre associated polynomials, like suggested in 
%   the reference below. This differs from the implementation in  Matlab's
%   'legendre' function (Version 2018b).
%
%   Parametrization of points on sphere:
%           Omega = [cos(phi) * sqrt(1 - mu^2);
%                    sin(phi) * sqrt(1 - mu^2);
%                    mu] 
% 
%   Implementation based on 
%       Evaluation of rotation matrices in the basis of real spherical
%       harmonics; Miguel A. Blanco, M. Florez, M. Bermejo; 1997; Journal
%       of Molecular Structure: THEOCHEM; Volume 419; 
%       doi: https://doi.org/10.1016/S0166-1280(97)00185-1
%
%   For a demonstration on how to use this function,
%   see also EVALONB
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

spatialDimension = 3;
nMoments = getNumberOfBasisFunctions(lmax, spatialDimension);

mu = pointsOnSphereMuPhi(1, :);
if (size(pointsOnSphereMuPhi, 1) == 1) || all(isnan(pointsOnSphereMuPhi(2, :)))
    phi = zeros(size(mu));
else
    phi = pointsOnSphereMuPhi(2, :);
end

P = evalLegendreAssociatedPolynomials(pointsOnSphereMuPhi, lmax);

C = cell(lmax, 1);
S = C;
[C{:}] = deal(cos(phi));
[S{:}] = deal(sin(phi));

for m = 2 : lmax
    C{m} = C{1} .* C{m-1} - S{1} .* S{m-1};
    S{m} = S{1} .* C{m-1} + C{1} .* S{m-1};
end

realSphericalHarmonicsAtPoints = zeros(nMoments, size(mu, 2));
N = @(l, m) sqrt( (2 * l + 1) / 4 / pi * factorial(l - m) / factorial(l + m) );

linIdx = 0 : nMoments - 1;
[l, m] = linearIdx2DegOrder(linIdx, spatialDimension);
for i = 1 : length(linIdx)
    if m(i) < 0
        realSphericalHarmonicsAtPoints(linIdx(i) + 1, :) = ...
            sqrt(2) * S{abs(m(i))} * N(l(i), abs(m(i))) .* P{l(i) + 1, abs(m(i)) + 1};
    end
    if m(i) == 0
        realSphericalHarmonicsAtPoints(linIdx(i) + 1, :) = ...
            N(l(i), m(i)) .* P{l(i) + 1, m(i) + 1};
    end
    if m(i) > 0
        realSphericalHarmonicsAtPoints(linIdx(i) + 1, :) = ...
            sqrt(2) .* C{m(i)} .* N(l(i), m(i)) .* P{l(i) + 1, m(i) + 1};
    end
end
end