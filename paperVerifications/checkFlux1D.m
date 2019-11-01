% CHECKFLUX1D Check the explicit representation of the flux matrix in 1D.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

clear
addpath(genpath('../PN2nd/'));
f = @(l)sqrt(1/(4 * l^2 + 8 * l + 3)) * (l + 1);
[~, ~, Tz] = fluxPNSphHarm(30, 1);

for l = 0 : 29
    fprintf('Should be zero: %d\n', f(l) - Tz(l + 1, l + 2));
end

