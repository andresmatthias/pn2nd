function kernelMatrix = assembleKernelMatrix2D(kernel, vx, vy, vz)
% ASSEMBLEKERNELMATRIX2D Assemble the kernel matrix for a 2D test case, 
%   where the entry i,j corresponds to an evaluation of the kernel function
%   at the directions Omega_i, Omega_j (Omega = [vx, vy, vz]).
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%
    n = length(vx);
    kernelMatrix = zeros(n, n);
    for j = 1:n
        kernelMatrix(:, j) = kernel.fun(vx, vy, vz, repmat(vx(j), 1, n),...
            repmat(vy(j), 1, n), repmat(vz(j), 1, n));
    end
end