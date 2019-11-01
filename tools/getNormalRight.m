function n = getNormalRight(edgeDirection)
% GETNORMALRIGHT Get the normal of a 2D edge (pointing to the right when
%   looking into the direction of the edge.)
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%
    n = [edgeDirection(2, :); -edgeDirection(1, :)];
    n = n ./ sqrt(sum(n.^2, 1));
end