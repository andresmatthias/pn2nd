function n = getNormalRight(edgeDirection)
% GETNORMALRIGHT Get the normal of a 2D edge (pointing to the right when
%   looking into the direction of the edge.)
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%
    n = [edgeDirection(2, :); -edgeDirection(1, :)];
    n = n ./ sqrt(sum(n.^2, 1));
end