function [indexSetTransportLhs, indexSetCollisionLhs, indexSetTransportRhs] = ...
    assembleRowDO2D(idx, par, mesh)
% ASSEMBLEROWDO2D Assemble rows of the discrete ordinates system matrix.
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
    assembleRowTransportTermDO2D(tupidx, par, mesh);

%% collision
indexSetCollisionLhs = assembleRowCollisionTermDO2D(tupidx, par, mesh);

end