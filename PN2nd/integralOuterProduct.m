function [I] = integralOuterProduct(fAtQuadrature, gAtQuadrature, weights)
% INTEGRALOUTERPRODUCT Compute the integral over the outer product f*g^T
%   for functions f in R^n and g in R^m.
% 
%   For a demonstration on how to use this function,
%   see also BASISFUNCTIONSTEST, FLUXMATRIXTEST, MISCTEST
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

if (size(fAtQuadrature, 2) ~= size(gAtQuadrature, 2)) || (size(fAtQuadrature, 2) ~= length(weights))
   error('Wrong input dimensions!') 
end
nf = size(fAtQuadrature, 1);
ng = size(gAtQuadrature, 1);

fAux = repmat(permute(fAtQuadrature, [1, 3, 2]), 1, ng, 1);
gAux = repmat(permute(gAtQuadrature, [3, 1, 2]), nf, 1, 1);
wAux = repmat(permute(weights, [1, 3, 2]), nf, ng, 1);
I = sum(fAux .* gAux .* wAux, 3);
I = round(I, 14);
end