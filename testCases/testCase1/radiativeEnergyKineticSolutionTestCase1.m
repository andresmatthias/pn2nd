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
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
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