function [tupidx] = lin2tupDO(linidx, nr_2ndIndex)
% LIN2TUPDO Convert linear index to tuple index.
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

    linidx = linidx(:)' - 1;
    tupidx = [floor(linidx / nr_2ndIndex) + 1;
              mod(linidx, nr_2ndIndex) + 1];
end