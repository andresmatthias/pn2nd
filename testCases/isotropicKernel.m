function [kernel] = isotropicKernel()
% ISOTROPICKERNEL Define the kernel function for isotropic scattering.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

kernel.fun = @(vx, vy, vz, wx, wy, wz) 1 / 4 / pi * ones(size(vx .* wx));
kernel.name = 'isotropic';
end

