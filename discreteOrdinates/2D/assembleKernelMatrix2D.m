function kernelMatrix = assembleKernelMatrix2D(kernel, vx, vy, vz)
% ASSEMBLEKERNELMATRIX2D Assemble the kernel matrix for a 2D test case, 
%   where the entry i,j corresponds to an evaluation of the kernel function
%   at the directions Omega_i, Omega_j (Omega = [vx, vy, vz]).
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%
    n = length(vx);
    kernelMatrix = zeros(n, n);
    for j = 1:n
        kernelMatrix(:, j) = kernel.fun(vx, vy, vz, repmat(vx(j), 1, n),...
            repmat(vy(j), 1, n), repmat(vz(j), 1, n));
    end
end