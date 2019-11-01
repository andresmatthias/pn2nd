function [linidx] = tup2linDO(tupidx, nr_2ndindex)
% TUP2LINDO Convert tuple index to linear index.
%   The solution has one degree of freedom for each cell / discrete
%   ordinate pair. 
%   linear index: label all degrees of freedom by numbers 1, 2, ...
%   tuple index: label all degrees of freedom as tuple of cell number and
%       discrete ordinate number.
%   nr_2ndIndex: number of discrete ordinates (second tuple entry)
%
%   For a demonstration on how to use this function,
%   see also MISCDISCRETEORDINATESTEST
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

if any(tupidx(2, :) > nr_2ndindex)
    error('something went wrong: ordinate index too large')
end
linidx = (tupidx(1, :) - 1 ) * nr_2ndindex + tupidx(2, :);
end