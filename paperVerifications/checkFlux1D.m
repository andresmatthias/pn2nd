% CHECKFLUX1D Check the explicit representation of the flux matrix in 1D.
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

clear
addpath(genpath('../PN2nd/'));
f = @(l)sqrt(1/(4 * l^2 + 8 * l + 3)) * (l + 1);
[~, ~, Tz] = fluxPNSphHarm(30, 1);

for l = 0 : 29
    fprintf('Should be zero: %d\n', f(l) - Tz(l + 1, l + 2));
end

