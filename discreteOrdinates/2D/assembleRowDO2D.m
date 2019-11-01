function [indexSetTransportLhs, indexSetCollisionLhs, indexSetTransportRhs] = ...
    assembleRowDO2D(idx, par, mesh)
% ASSEMBLEROWDO2D Assemble rows of the discrete ordinates system matrix.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

tupidx = lin2tupDO(idx, par.nDiscreteOrdinates); % 1: physical (element), 2: directional (ordinate)

%% transport
[indexSetTransportLhs, indexSetTransportRhs] = ...
    assembleRowTransportTermDO2D(tupidx, par, mesh);

%% collision
indexSetCollisionLhs = assembleRowCollisionTermDO2D(tupidx, par, mesh);

end