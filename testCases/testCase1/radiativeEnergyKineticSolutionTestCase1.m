function radiativeEnergy = radiativeEnergyKineticSolutionTestCase1(...
    externalSourceZMin, externalSourceZMax, zGrid, sigma_aConst)
% RADIATIVEENERGYKINETICSOLUTIONTESTCASE1 Evaluate the analytic solution of 
%   the kinetic problem for test case 1 at given grid points and boundary 
%   sources.
%
%   domain: [0, 1]
%   sigma_a independent of z
%   vz * d_z I = - sigma_a * I
%   I(z=0, vx, vy, vz) = S_left(vx, vy, vz)
%   I(z=1, vx, vy, vz) = S_right(vx, vy, vz) 
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%
I = @(z, vx, vy, vz) kineticSolution(z, vx, vy, vz);

maxExactDegree = 50;
[weights, mu, phi] = sphericalQuadratureFull(maxExactDegree);
Omega = sphereParam2Cartesian(mu, phi);
radiativeEnergy = zeros(length(zGrid), 1);
for k = 1 : length(zGrid)
   radiativeEnergy(k) = sum(weights .* I(zGrid(k), Omega(1, :), Omega(2, :), Omega(3, :))); 
end

function I = kineticSolution(z, vx, vy, vz)
    I = zeros(size(vz));
    idx = vz > 1e-14;
    I(idx) = exp(- sigma_aConst ./ vz(idx) * z) .* externalSourceZMin(vx(idx), vy(idx), vz(idx));
    
    idx = vz < -1e-14;
    I(idx) = exp(sigma_aConst ./ vz(idx) * (1 - z)) .* externalSourceZMax(vx(idx), vy(idx), vz(idx));
end
end