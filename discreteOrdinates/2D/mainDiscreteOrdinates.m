function [PHI, I, par, A_upwind, rhs_upwind] = mainDiscreteOrdinates(par, mesh, maxExactDegree)
% MAINDISCRETEORDINATES2D Compute the discrete ordinates solution of a
%   2D test case defined by its parameter set.
%   Result is interpreted as solution of a problem symmetric along the 
%   z-axis (i.e., as a 3D solution with dz= 0).
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

    %% attenuation coefficient
    par.sigma_t = @(x) par.sigma_a(x) + par.sigma_s(x);
    %% check mesh
    boundaryEdgesByPointIdx = mesh.edgeList(:, mesh.boundaryEdgesIdx2Id.edgeIdx);
    sortedBoundaryEdges = sort(boundaryEdgesByPointIdx, 1);
    if norm(boundaryEdgesByPointIdx - sortedBoundaryEdges) > 0
        error('Check mesh properties: indices of boundary edge need to be sorted')
    end
    [~, changeFlag] = orderConnectivityListCcw(mesh.points, mesh.connectivityList);
    if changeFlag
       error('Connectivity list not ordered counter clockwise: abort') 
    end
    %% generate discrete ordinates
    par.ordinates = getOrdinates(maxExactDegree, par.spatialDimension, par.kernel.name);
    par.nDiscreteOrdinates = length(par.ordinates.weights);
    %% assemble kernel matrix
    par.kernelMatrix = assembleKernelMatrix2D(par.kernel,...
        par.ordinates.Cartesian(1, :), par.ordinates.Cartesian(2, :), par.ordinates.Cartesian(3, :));
    
    %% assemble the system matrix and right hand side
    fprintf('\nAssemble rows \n')
    % indices grouped in blocks, where each block corresponds to a single
    % physical element's index, and contains all directional indices at the specified
    % physical point, i.e., 1:(t_1,s_1),2:(t_1,s_2),...,No:(t_1,s_No),
    % No+1:(t_2,s_1), ... 2*No:(t_2,s_No),...
    
    indices = 1 : mesh.nElements * par.nDiscreteOrdinates;
    %% capacity check
    fprintf('\n\n no degrees of freedom: %d \n no discrete ordinates: %d \n\n',...
        length(indices), par.nDiscreteOrdinates)
    %%
    fAux = @(idx) assembleRowDO2D(idx, par, mesh);
    
    [indexSetTransportLhs, indexSetCollisionLhs, indexSetRhs] =...
            fAux(indices);
    
    indexSetLhs = cat(1, indexSetTransportLhs, indexSetCollisionLhs);

    clear index_set_transport index_set_collision 
    
    %%
    fprintf('Create sparse matrix')
    A_upwind = sparse(indexSetLhs(:, 1), indexSetLhs(:, 2), indexSetLhs(:, 3), ...
        mesh.nElements * par.nDiscreteOrdinates, ...
        mesh.nElements * par.nDiscreteOrdinates); % sums up double entries
    
    fprintf('Create sparse rhs')
    rhs_upwind = sparse(indexSetRhs(:, 1), indexSetRhs(:, 2), ...
        indexSetRhs(:, 3), mesh.nElements * par.nDiscreteOrdinates, 1);
    
    fprintf('\n Solve linear system \n')
    I = A_upwind \ rhs_upwind;    
     
    fprintf('\n Integrate to obtain density \n')
    
    %% quadrature to obtain phi from I 
    PHI = zeros(size(mesh.connectivityList, 2), 1);
    for k = 1:length(PHI)
        idx = tup2linDO([repmat(k, 1, par.nDiscreteOrdinates); ...
            (1 : par.nDiscreteOrdinates)], par.nDiscreteOrdinates);
        PHI(k) = (par.ordinates.weights *  I(idx));
    end
    fprintf('\n finished :) \n')
end
