function kernelMatrix = assembleKernelMatrix1D(kernel, mu, phi)
% ASSEMBLEKERNELMATRIX1D Assemble the kernel matrix for a 1D test case, 
%   where the entry i,j corresponds to an evaluation of the kernel function
%   at the directions mu_i, mu_j .
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%
    n = length(mu);
    kernelMatrix = zeros(n, n);
    Omega = sphereParam2Cartesian(mu, phi);
    for j = 1 : n
        vx = repmat(Omega(1, j), 1, n);
        vy = repmat(Omega(2, j), 1, n);
        vz = repmat(Omega(3, j), 1, n);
        wx = Omega(1, :);
        wy = Omega(2, :);
        wz = Omega(3, :);
        kernelMatrix(:, j) = kernel.fun(wx, wy, wz, vx, vy, vz);
    end
end