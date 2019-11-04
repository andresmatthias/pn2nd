function [kernel] = isotropicKernel()
% ISOTROPICKERNEL Define the kernel function for isotropic scattering.
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

kernel.fun = @(vx, vy, vz, wx, wy, wz) 1 / 4 / pi * ones(size(vx .* wx));
kernel.name = 'isotropic';
end

