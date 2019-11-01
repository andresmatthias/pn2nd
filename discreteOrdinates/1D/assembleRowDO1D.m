function [indexSetTransportLhs, indexSetCollisionLhs, indexSetTransportRhs] = ...
    assembleRowDO1D(idx, par)
% ASSEMBLEROWDO1D Assemble rows of the discrete ordinates system matrix.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

tupidx = lin2tupDO(idx, par.nDiscreteOrdinates); % 1: physical (element), 2: directional (ordinate)

%% transport
[indexSetTransportLhs, indexSetTransportRhs] = ...
    assembleRowTransportTermDO1D(tupidx, par);

%% collision
indexSetCollisionLhs = assembleRowCollisionTermDO1D(tupidx, par);

end