function radiativeEnergy = radiativeEnergyKineticSolutionTestCase3(zGrid)
% RADIATIVEENERGYKINETICSOLUTIONTESTCASE3 Evaluate the analytic solution of 
%   the kinetic problem for test case 3 at given grid points and boundary 
%   sources.
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%
I = @(z, vx, vy, vz) vz - z .* (z + 2) / 3 + 2;
dzI = @(z, vx, vy, vz) - 2 / 3 * z - 2 / 3;
sigma_s = @(z) 1 + z;
sigma_a = @(z) 0;
radiativeEnergy = 2 * pi * (4 - 2 / 3 * zGrid .* (zGrid + 2));

kernel = @(vx, vy, vz, wx, wy, wz) 1 / 8 / pi * ((vz - 1).*(wz - 1) + (vz + 1).*(wz + 1));
maxExactDegree = 13;
[weights, mu, phi] = sphericalQuadratureFull(maxExactDegree);
Omega = sphereParam2Cartesian(mu, phi);

%% check solution
% mu * dz I = -(sigma_a + sigma_s) * I + sigma_s * int_{4*pi} kernel(mu, mu') * I(mu') dmu'
checkFlag = 0;
for j = 1 : length(mu)
    for k = 1 : length(zGrid)
        tmp = mu(j) .* dzI(zGrid(k), Omega(1, j), Omega(2, j), Omega(3, j))...
            + (sigma_a(zGrid(k)) + sigma_s(zGrid(k))) * I(zGrid(k), Omega(1, j), Omega(2, j), Omega(3, j)) ...
            - sigma_s(zGrid(k)) * sum(weights .* kernel(Omega(1, j), Omega(2, j), Omega(3, j), Omega(1, :), Omega(2, :), Omega(3, :))...
             .* I(zGrid(k), Omega(1, :), Omega(2, :), Omega(3, :)));
        if abs(tmp) > 1e-13
           checkFlag = checkFlag + 1; 
        end
    end
end
if checkFlag > 0
    fprintf('Something wrong with analytical kinetic solution in test case 1D 3.\n')
else
    fprintf('Analytical kinetic solution in test case 1D 3 correct.\n')
end
end