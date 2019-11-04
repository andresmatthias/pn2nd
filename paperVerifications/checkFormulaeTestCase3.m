% CHECKFORMULASTESTCASE3 Check formulas regarding test case 3.
%
%   rotational symmetry around z-Axis:
%   mu * d_z I = (1 + z)int_{4pi} k(Omega,Omega') I(Omega') d Omega
%   with
%   k(Omega, Omega') = 1 / 8 / pi * ((mu-1)*(mu'-1) + (mu+1)*(mu'+1))
%   I(z=0) = mu + 2 (left boundary, zero reflection)
%   I(z=1) = mu + 1 (right boundary, zero reflection)
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%


clear
addpath(genpath('../PN2nd/'));
addpath('../testCases/testCase3/')
%% define candidate for the solution of the kinetic equation
syms mu z phi mup
I = mu - z * (z + 2) / 3 + 2;
Ip = subs(I, mu, mup);
k = sym(1 / 8) / pi * ((mu - 1) * (mup - 1) + (mu + 1) * (mup + 1));
collisionIntegral = int(int(symfun(k * Ip, [mup, phi]), mup, -1, 1), phi, 0, 2*pi);

fprintf('\nKinetic problem\n---------------\n')
fprintf('Equation on domain , should be zero:\t\t %s\n', simplify(diff(I, z) * mu + (1+z) * I - (1 + z) * collisionIntegral));

fprintf('Left boundary (mu >= 0), should be mu + 2:\t %s\n', subs(I, z, 0))
fprintf('Right boundary (mu < 0), should be mu + 1:\t %s\n', subs(I, z, 1))

fprintf('I:\n')
pretty(I)

