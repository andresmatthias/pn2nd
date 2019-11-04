function [] = plotPwConstantFunction1D(fPwConst, grid, varargin)
% PLOTPWCONSTANTFUNCTION1D Plot a piecewise constant function in 1D.
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

if (length(grid)-1) ~= length(fPwConst)
   error('sth wrong with dimensions') 
end
for k = 1:length(grid)-1
    plot(grid(k:k+1),[fPwConst(k), fPwConst(k)], varargin{:})
    hold on
end
end