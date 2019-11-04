function [kernel] = bilinearKernel()
% BILINEARKERNEL Define the kernel function for bilinear scattering.
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%
kernel.fun = @(vx, vy, vz, wx, wy, wz) 1 / 8 / pi * ((vz - 1) .* (wz - 1) + (vz + 1) .* (wz + 1));
kernel.name = 'bilinear';
end

