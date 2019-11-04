function idxRed = getReducedIdx(l, m, spatialDimension)
% GETREDUCEDIDX Return indices of those basis functions (real spherical
%   harmoncis) included after reduction due to symmetry (spatialDimension).
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

switch spatialDimension 
    case 1
        % reduction as no dependence on phi
        idxRed = (m == 0);
    case 2
        % reduction as no dependence on mu
        idxRed = (mod(l + m, 2) == 0);
    case 3
        idxRed = true(size(l));  % take care of Matlab shift
end
end