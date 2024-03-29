function [OmegaRef] = reflectionAtSurface(outerNormalVector, Omega)
% REFLECTIONATSURFACE Compute the reflected vector at a surface in 3D
%   (given by its outer normal vector).
% 
%   For a demonstration on how to use this function,
%   see also MISCTEST
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

if (size(outerNormalVector, 1) ~= 3) || (size(Omega, 1) ~= 3)
   error('Wrong input dimensions!') 
end
scp = sum(Omega .* outerNormalVector, 1);
if any(scp > 1e-15)
   error('Omega does not point inside domain, w.r.t. outer normal vector!') 
end
OmegaRef = Omega - 2 * repmat(scp, 3, 1) .* outerNormalVector; 
end