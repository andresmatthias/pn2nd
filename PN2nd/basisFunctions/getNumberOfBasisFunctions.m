function nMoments = getNumberOfBasisFunctions(N, spatialDimension)
% GETNUMBEROFBASISFUNCTIONS  Get number of reduced angular basis functions
%   up to degree N, depending on the spatial dimension.
%
%   For a demonstration on how to use this function,
%   see also BASISFUNCTIONSTEST
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

switch spatialDimension
    case 1
        nMoments = N + 1;
    case 2
        nMoments = N^2 / 2 + 3 / 2 * N + 1;
    case 3
        nMoments = N^2 + 2 * N + 1;
    otherwise
        error('Invalid spatial dimension!') 
end
end