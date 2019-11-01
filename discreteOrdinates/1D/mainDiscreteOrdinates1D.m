function [radiativeEnergy] = mainDiscreteOrdinates1D(par)
% MAINDISCRETEORDINATES1D Compute the discrete ordinates solution of a
%   1D test case defined by its parameter set.
%   Result is interpreted as solution of 1D slab geometry (i.e., as a
%   3D solution with dx = dy = 0).
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

%% attenuation coefficient
par.sigma_t = @(x) par.sigma_a(x) + par.sigma_s(x);
%% generate discrete ordinates, reduced
par.ordinates = getOrdinates(par.maxExactDegree, par.spatialDimension, par.kernel.name);
par.nDiscreteOrdinates = length(par.ordinates.mu);
fprintf('(Reduced) number of discrete ordinates: %d \n', par.nDiscreteOrdinates);
%% compute element areas
par.areaElement = diff(par.grid);
par.nElements = par.nGridPoints - 1;
%% assemble kernel matrix
par.kernelMatrix = assembleKernelMatrix1D(par.kernel, par.ordinates.mu, par.ordinates.phi);
%% assemble the system matrix and right hand side
fprintf('\nAssemble rows \n')
% indices grouped in blocks, where each block corresponds to a single
% physical element's index, and contains all directional indices at the specified
% physical point, i.e., 1:(t_1,s_1),2:(t_1,s_2),...,No:(t_1,s_No),
% No+1:(t_2,s_1), ... 2*No:(t_2,s_No),...

par.idxMax = par.nElements * par.nDiscreteOrdinates;
indices = 1 : par.idxMax;
%%
faux = @(idx) assembleRowDO1D(idx, par);

[indexSetTransportLhs, indexSetCollisionLhs, indexSetRhs] = faux(indices);

indexSetLhs = cat(1, indexSetTransportLhs, indexSetCollisionLhs);

clear index_set_transport index_set_collision 

%%
fprintf('Create sparse matrix')
lhsUpwind = sparse(indexSetLhs(:, 1), indexSetLhs(:, 2), indexSetLhs(:, 3), ...
    par.idxMax, par.idxMax); % sums up double entries

fprintf('Create sparse rhs')
rhsUpwind = sparse(indexSetRhs(:, 1), indexSetRhs(:, 2), ...
    indexSetRhs(:, 3), par.idxMax, 1);

fprintf('\n Solve linear system \n')
I = lhsUpwind \ rhsUpwind;  

fprintf('\n Integrate to obtain density \n')

%% quadrature to obtain radiative energy (3D) from intensity
radiativeEnergy = zeros(1, par.nElements);
for k = 1:length(radiativeEnergy)
    idx = tup2linDO([repmat(k, 1, par.nDiscreteOrdinates); ...
        (1 : par.nDiscreteOrdinates)], par.nDiscreteOrdinates);
    radiativeEnergy(k) = par.ordinates.weights *  I(idx);
end
fprintf('\n finished :) \n')
end