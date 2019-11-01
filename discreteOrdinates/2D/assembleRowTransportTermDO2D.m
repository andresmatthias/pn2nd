function [indexSetLhs, indexSetRhs] = assembleRowTransportTermDO2D(tupidx, par, mesh)
% ASSEMBLEROWTRANSPORTTERMDO2D Assemble rows of the discrete ordinates
%   system matrix regarding the transport term in the original kinetic
%   equation.
% 
%   tupidx: first row: element indices
%           second row: ordinate indices
% 
%#$&replaceMe&$#%

% for transfer part: entry in (i,j) can only be coupled with
% (i,j),(i1,j),(i2,j),(i3,j), where i1,i2,i3 are indices of
% neighbouring elements to i; and, if i is a boundary element, due to
% reflection with (i,j'), where  j' corresponds to the reflected direction
% for reflection: per reflected ordinate we get three indices due to
% possible interpolation

%% handles the following part:
%
% $$\frac{1}{|T_i|}\int_{T_i} \Omega_j\cdot \nabla I_j^i$$
%
nIdx = size(tupidx, 2);

areaT = mesh.elementVolume(tupidx(1, :));

% neighbouring elements (= NaN if boundary edge)
% ordered opposite of points in connectivity_list
idx_T_adjacent = mesh.elementNeighbors(:, tupidx(1, :));

% NaN if exterior edge
% 1: (el, ord), 2: (el_a, ord), 3: (el_b, ord), 4:(el_c, ord),
% 5: (el, ord_a'), 6: (el, ord_b'), 7: (el, ord_c')
% indices of reflected indices will be assigned below
adjacencyThroughEdgeBC = [idx_T_adjacent(1, :);
                             tupidx(2, :)];
adjacencyThroughEdgeCA = [idx_T_adjacent(2, :);
                             tupidx(2, :)];
adjacencyThroughEdgeAB = [idx_T_adjacent(3, :);
                             tupidx(2, :)];

adjacentIndices = [tup2linDO(tupidx, par.nDiscreteOrdinates);
    tup2linDO(adjacencyThroughEdgeBC, par.nDiscreteOrdinates);
    tup2linDO(adjacencyThroughEdgeCA, par.nDiscreteOrdinates);
    tup2linDO(adjacencyThroughEdgeAB, par.nDiscreteOrdinates);
    NaN(9, nIdx)]; % 9 = 3 x 3; 
% one entry for each possible reflection at edge, and three for interpolation
% of each

matrixEntries = zeros(13, nIdx);
constantPart = zeros(1, nIdx); % results from boundary terms

Omega = par.ordinates.Cartesian(:, tupidx(2, :));
%% run over edges of element
% ordered for each element according to neighbors, i.e. opposite of points
% in connectivity list
for k = 1:3
    edgeLen_k = mesh.edgeLengths(mesh.edgesPerElement(k, tupidx(1, :)));
    outerNormal_k = squeeze(mesh.edgeNormalsPerElement(:, k, tupidx(1, :)));
    
    scalarProduct_k  = dot(Omega(1:2, :), outerNormal_k, 1);
    factor = edgeLen_k ./ areaT .* scalarProduct_k;
    
    %% flow outward, take radiation in current cell, flowing outward
    idxOutwards_k = (scalarProduct_k>=0);
    matrixEntries(1, idxOutwards_k) = matrixEntries(1, idxOutwards_k) + factor(idxOutwards_k);
    
    %%  flow inward
    idxInwards_k = (scalarProduct_k<0);
    
    % interior edge;
    idxInwardsInterior_k = (idxInwards_k & ~isnan(idx_T_adjacent(k, :)));
    matrixEntries(k + 1, idxInwardsInterior_k) = factor(idxInwardsInterior_k);
    
    % boundary edges
    idxInwardsBoundary_k = (idxInwards_k & isnan(idx_T_adjacent(k, :)));
    inwardsBoundaryIDs_k = mesh.boundaryIdsPerElement(k, tupidx(1, idxInwardsBoundary_k));
    rho = par.reflectivity(inwardsBoundaryIDs_k);
    
    aux = unique(inwardsBoundaryIDs_k);
    SExt = zeros(size(inwardsBoundaryIDs_k));
    for j = aux
       idxAux = inwardsBoundaryIDs_k == j;
       OmegaAux = Omega(:, idxInwardsBoundary_k);
       SExt(idxAux) = par.externalSource{j}(OmegaAux(1, idxAux),...
           OmegaAux(2, idxAux), OmegaAux(3, idxAux));
    end
    
    % external inflow
    % add, as same  tuple index could occur for two different boundary
    % edges
    constantPart(idxInwardsBoundary_k) = constantPart(idxInwardsBoundary_k)...
        + factor(idxInwardsBoundary_k) .* (1 - rho) .* SExt;
    
    % reflection
    [outgoinOrdinateIdx, reflectionInterpWeights] = getOutgoingOrdinateIdxInterp2D(...
        tupidx(2, idxInwardsBoundary_k), outerNormal_k(:,...
        idxInwardsBoundary_k), par.ordinates);
    
    adjacencyThroughEdge_k = [repmat(tupidx(1, idxInwardsBoundary_k), 1, 3);
                                outgoinOrdinateIdx(1, :),...
                                outgoinOrdinateIdx(2, :),...
                                outgoinOrdinateIdx(3, :)];
    tmp = tup2linDO(adjacencyThroughEdge_k, par.nDiscreteOrdinates);
    normalIn = sum(idxInwardsBoundary_k);
    tmp = [tmp(1 : normalIn); tmp(normalIn + 1 : 2 * normalIn); tmp(2 * normalIn + 1 : end)];
    adjacentIndices((k - 1) * 3 + 4 + 1 : k * 3 + 4, idxInwardsBoundary_k) = tmp;
    matrixEntries((k - 1) * 3 + 4 + 1 : k * 3 + 4, idxInwardsBoundary_k) = ...
        repmat(factor(idxInwardsBoundary_k) .* rho, 3, 1) .* reflectionInterpWeights;
    
end

%% create index list for sparse representation
% throw out indices corresponding to exterior edges / ghost cells
% instead of throwing out: add zero
indexSetLhs = [reshape(repmat(1 : nIdx, 13, 1), [], 1), adjacentIndices(:), matrixEntries(:)];
idxLhs = ~isnan(indexSetLhs(:, 2));
indexSetLhs = indexSetLhs(idxLhs, :);
     
allLinIdx = 1 : nIdx;
idxRhs = (constantPart ~= 0);
% attention to minus! 
indexSetRhs = [allLinIdx(idxRhs)', ones(sum(idxRhs), 1), -constantPart(idxRhs)'];

end

