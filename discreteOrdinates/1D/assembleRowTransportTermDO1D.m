function [indexSetLhs, indexSetRhs] = assembleRowTransportTermDO1D(tupidx, par)
% ASSEMBLEROWTRANSPORTTERMDO1D Assemble rows of the discrete ordinates
%   system matrix regarding the transport term in the original kinetic
%   equation.
% 
%   tupidx: first row: element indices
%           second row: ordinate indices
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

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

areaT = par.areaElement(tupidx(1, :));

% - same cell and ordinate
% - left cell and same ordinate
% - right cell and same ordinate
% - same cell but reflected ordinate
adjacentIndices = [tup2linDO(tupidx, par.nDiscreteOrdinates);
            tup2linDO([tupidx(1, :) - 1; tupidx(2, :)], par.nDiscreteOrdinates);
            tup2linDO([tupidx(1, :) + 1; tupidx(2, :)], par.nDiscreteOrdinates);
            NaN(1, nIdx)]; 
adjacentIndices(adjacentIndices < 1 ) = nan;      
adjacentIndices(adjacentIndices >  par.idxMax) = nan;      

% one entry for each possible reflection at edge, and 3 for interpolation
% of each

matrixEntries = zeros(4, nIdx);
constantPart = zeros(1, nIdx); % results from boundary terms

mu = par.ordinates.mu(tupidx(2, :));

%% interior
knudsenNr = 1;
matrixEntries(1, :) = knudsenNr ./ areaT .* abs(mu);
matrixEntries(2, :) = - knudsenNr ./ areaT .* max(mu, 0);
matrixEntries(3, :) = knudsenNr ./ areaT .* min(mu, 0);

%% boundary
leftBoundaryIdx = (tupidx(1, :) == 1) & (mu > 0);
rightBoundaryIdx = (tupidx(1, :) == par.nElements) & (mu < 0);

% [outgoingOrdinateLeftIdx] = getOutgoingOrdinateIdx1D(...
%     tupidx(2, leftBoundaryIdx), par.ordinates.mu);
% [outgoingOrdinateRightIdx] = getOutgoingOrdinateIdx1D(...
%     tupidx(2, rightBoundaryIdx), par.ordinates.mu);

[outgoingOrdinateLeftIdx] = getOutgoingOrdinateIdx1D(...
    tupidx(2, leftBoundaryIdx), par.ordinates);
[outgoingOrdinateRightIdx] = getOutgoingOrdinateIdx1D(...
    tupidx(2, rightBoundaryIdx), par.ordinates);

adjacentIndices(4, leftBoundaryIdx) = tup2linDO([ones(1, sum(leftBoundaryIdx));...
                                                  outgoingOrdinateLeftIdx], par.nDiscreteOrdinates);
adjacentIndices(4, rightBoundaryIdx) = tup2linDO([par.nElements * ones(1, sum(rightBoundaryIdx));...
                                                  outgoingOrdinateRightIdx], par.nDiscreteOrdinates);

% max / min in the following redundant, but for better readability we keep
% it
% 1: left
% 2: right
matrixEntries(4, leftBoundaryIdx) = - knudsenNr ./ areaT(leftBoundaryIdx) ...
    * par.reflectivity(1) .* max(mu(leftBoundaryIdx), 0);
matrixEntries(4, rightBoundaryIdx) = knudsenNr ./ areaT(rightBoundaryIdx) ...
    * par.reflectivity(2) .* min(mu(rightBoundaryIdx), 0);

%% boundary external source
% ID = 1: left
% ID = 2: right
sourceLeft = par.externalSource{1};
vxyDummy = zeros(1, sum(leftBoundaryIdx));
constantPart(leftBoundaryIdx) = -knudsenNr ./ areaT(leftBoundaryIdx) ...
    .* (1 - par.reflectivity(1)) .* sourceLeft(vxyDummy, vxyDummy, mu(leftBoundaryIdx))...
    .* max(mu(leftBoundaryIdx), 0);

sourceRight = par.externalSource{2};
vxyDummy = zeros(1, sum(rightBoundaryIdx));
constantPart(rightBoundaryIdx) = knudsenNr ./ areaT(rightBoundaryIdx) ...
    .* (1 - par.reflectivity(2)) .* sourceRight(vxyDummy, vxyDummy, mu(rightBoundaryIdx))...
    .* min(mu(rightBoundaryIdx), 0);
                                              
%% create index list for sparse representation
% throw out indices corresponding to exterior edges / ghost cells
% instead of throwing out: add zero
indexSetLhs = [reshape(repmat(1 : nIdx, 4, 1), [], 1), adjacentIndices(:), matrixEntries(:)];
idx_lhs = ~isnan(indexSetLhs(:, 2));
indexSetLhs = indexSetLhs(idx_lhs, :);
     
allLinIdx = 1 : nIdx;
idxRhs = (constantPart ~= 0);
% attention to minus! 
indexSetRhs = [allLinIdx(idxRhs)', ones(sum(idxRhs), 1), -constantPart(idxRhs)'];

end

