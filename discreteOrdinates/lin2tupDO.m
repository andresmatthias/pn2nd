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
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

    linidx = linidx(:)' - 1;
    tupidx = [floor(linidx / nr_2ndIndex) + 1;
              mod(linidx, nr_2ndIndex) + 1];
end