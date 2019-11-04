function idxOdd = linearIdxOfOddBasis(N, spatialDimension)
% LINEARIDXOFODDBASIS Get linear indices of all basis functions (real 
%   spherical harmonics) with odd degree <= N.
%
%   Depending on spatial dimension obtain the odd functions out of the
%   reduced set. 
%
%   For a demonstration on how to use this function,
%   see also BASISFUNCTIONSTEST
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

nMoments = getNumberOfBasisFunctions(N, spatialDimension);
linIdx = 0 : nMoments - 1;
[l, ~] = linearIdx2DegOrder(linIdx, spatialDimension);
idxOdd = (mod(l, 2) == 1);
idxOdd = find(idxOdd) - 1;
end