function [kernel] = bilinearKernel()
% BILINEARKERNEL Define the kernel function for bilinear scattering.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%
kernel.fun = @(vx, vy, vz, wx, wy, wz) 1 / 8 / pi * ((vz - 1) .* (wz - 1) + (vz + 1) .* (wz + 1));
kernel.name = 'bilinear';
end

