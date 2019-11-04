function [indexSetTransportLhs, indexSetCollisionLhs, indexSetTransportRhs] = ...
    assembleRowDO1D(idx, par)
% ASSEMBLEROWDO1D Assemble rows of the discrete ordinates system matrix.
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

tupidx = lin2tupDO(idx, par.nDiscreteOrdinates); % 1: physical (element), 2: directional (ordinate)

%% transport
[indexSetTransportLhs, indexSetTransportRhs] = ...
    assembleRowTransportTermDO1D(tupidx, par);

%% collision
indexSetCollisionLhs = assembleRowCollisionTermDO1D(tupidx, par);

end