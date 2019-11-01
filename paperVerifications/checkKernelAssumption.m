function [] = checkKernelAssumption()
% CHECKKERNELASSUMPTION Check, that several well-known kernels (see 
%   http://www.eugenedeon.com/hitchhikers don't produce a drift term in the
%   second-order formulation of the PN equations.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

% isotropic scattering
k = @(vx, vy, vz, wx, wy, wz) 1 / 4 / pi * ones(size(vx));

%Eddington kernel
% b = rand(1);
% k = @(vx, vy, vz, wx, wy, wz) (1 + b * (vx .* wx + vy .* wy + vz .* wz)) / 4 / pi;

%Rayleigh kernel
% k = @(vx, vy, vz, wx, wy, wz) 3 / 16 / pi * (vx .* wx + vy .* wy + vz .* wz).^2;

%Kagiwada-Kalaba
% b = 0.5;
% k = @(vx, vy, vz, wx, wy, wz) b / 2 / pi ./ (1 - b * (vx .* wx + vy .* wy + vz .* wz)) ./ log((b + 1) / (1 - b));

%vonMisesFischer
% c = 15;
% k = @(vx, vy, vz, wx, wy, wz) c* exp(c * (vx .* wx + vy .* wy + vz .* wz)) * csch(c) / 4 / pi;

addpath(genpath('../PN2nd/'));
addpath('../discreteOrdinates/2D/');
[weights, mu, phi] = sphericalQuadratureFull(30);
Omega = sphereParam2Cartesian(mu, phi);


% Do some checks:
disp('----Symmetry----')
sym_flag = 0;
for i = 1 : size(mu, 2)
   for j = 1 : size(mu, 2)
        if  abs(k(Omega(1, i), Omega(2, i), Omega(3, i), Omega(1, j), Omega(2, j), Omega(3, j)) ...
                - k(Omega(1, j), Omega(2, j), Omega(3, j), Omega(1, i), Omega(2, i), Omega(3, i))) > 1e-13
        sym_flag = 1;
        disp('symmetry failed!') 
        break
        end
   end
end
if sym_flag == 0
   disp('symmetry fine!')
end

disp('-----Normalization-----')
norm_flag = 0;
for i = 1 : size(mu, 2)
    if abs(1 - sum(weights .* k(Omega(1, i), Omega(2, i), Omega(3, i), Omega(1, :), Omega(2, :), Omega(3, :)))) > 1e-13
       norm_flag = 1;
       disp('normalization failed!')
       break
    end
end
if norm_flag == 0
   disp('normalization fine!') 
end

disp('-----Positivity-----')
pos_flag = 0;
for i = 1 : size(mu, 2)
   for j = 1 : size(mu, 2)
       if k(Omega(1, i), Omega(2, i), Omega(3, i), Omega(1, j), Omega(2, j), Omega(3, j)) ...
                - k(Omega(1, j), Omega(2, j), Omega(3, j), Omega(1, i), Omega(2, i), Omega(3, i)) < 0
          pos_flag = 1;
          disp('positivity failed!')
          break
       end
   end
end
if pos_flag == 0
   disp('positivity fine!') 
end

% drift term:
spatialDimension = 3;
N = 11;
b = evalONB([mu; phi], N, spatialDimension);
idxEven = linearIdxOfEvenBasis(N, spatialDimension);
idxOdd = linearIdxOfOddBasis(N, spatialDimension);
be = b(idxEven + 1, :);
bo = b(idxOdd + 1, :);
Sigma_eo = zeros(length(idxEven), length(idxOdd));
kernel.fun = k;
kernelMatrix = assembleKernelMatrix2D(kernel, Omega(1, :), Omega(2, :), Omega(3, :));
for i = 1 : size(Sigma_eo, 1)
    for j = 1 : size(Sigma_eo, 2)
         Sigma_eo(i, j) = (be(i, :) .* weights) * kernelMatrix * (bo(j, :) .* weights)';
    end
end
disp('\nmax entry in Sigma_eo (should be (close to) zero)')
disp(max(abs(Sigma_eo(:))))
end


