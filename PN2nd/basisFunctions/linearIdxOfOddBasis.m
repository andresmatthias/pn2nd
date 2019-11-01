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
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

nMoments = getNumberOfBasisFunctions(N, spatialDimension);
linIdx = 0 : nMoments - 1;
[l, ~] = linearIdx2DegOrder(linIdx, spatialDimension);
idxOdd = (mod(l, 2) == 1);
idxOdd = find(idxOdd) - 1;
end