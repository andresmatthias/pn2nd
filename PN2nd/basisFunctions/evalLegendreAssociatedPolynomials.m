function P = evalLegendreAssociatedPolynomials(pointsOnSphereMuPhi, lmax)
% EVALLEGENDREASSOCIATEDPOLYNOMIALS Evaluate all Legendre associated
%   polynomials up to degree lmax on at given points on sphere.
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
%   The output cell is ordered as follows:
%       (l=0, m=0)
%       (l=1, m=0) | (l=1, m=1)
%       (l=2, m=0) | (l=2, m=1) | (l=2, m=2) 
%       (l=3, m=0) | ...
% 
%   Implementation based on 
%       Evaluation of rotation matrices in the basis of real spherical
%       harmonics; Miguel A. Blanco, M. Florez, M. Bermejo; 1997; Journal
%       of Molecular Structure: THEOCHEM; Volume 419; 
%       doi: https://doi.org/10.1016/S0166-1280(97)00185-1
%
%   For a demonstration on how to use this function,
%   see also BASISFUNCTIONSTEST
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

mu = pointsOnSphereMuPhi(1, :);

P = cell(lmax + 1);
l = 0;
P{l + 1, l + 1} = ones(size(mu));

for l = 1 : lmax
    P{l + 1, l + 1} = (2 * l - 1) * sqrt(1 - mu.^2) .* P{l - 1 + 1, l - 1 + 1};
    P{l + 1, l - 1 + 1} = (2 * (l - 1) + 1) * mu .* P{l - 1 + 1, l - 1 + 1};
end

for l = 2 : lmax
    for m = 0 : l - 2
        P{l + 1, m + 1} = ((2 * l - 1) * mu .* P{l - 1 + 1, m + 1}...
            - (l + m - 1) * P{l - 2 + 1, m + 1}) / (l - m);
    end
end

end