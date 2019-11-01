function [] = plotPwConstantFunction1D(fPwConst, grid, varargin)
% PLOTPWCONSTANTFUNCTION1D Plot a piecewise constant function in 1D.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

if (length(grid)-1) ~= length(fPwConst)
   error('sth wrong with dimensions') 
end
for k = 1:length(grid)-1
    plot(grid(k:k+1),[fPwConst(k), fPwConst(k)], varargin{:})
    hold on
end
end