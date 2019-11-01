function [l, m] = linearIdx2DegOrder(linIdx, spatialDimension)
% LINEARIDX2DEGORDER Convert the linear index (when serially numbering all
%   basis functions, depending on symmetry assumptions (spatialDimension))
%   to a tuple with the degree and order of the corresponding 
%   real spherical harmonic (S_l^m: l:degree, m: order).
%
%   3D (-l<=m<=l)     2D(l+m even)  1D(m=0)
%   n: (l, m)         n: (l, m)     n: (l, m)
%   0: (0, 0)         0: (0, 0)     0: (0, 0)
%   1: (1, -1)        1: (1, -1)    
%   2: (1, 0)                       1: (1, 0)
%   3: (1, 1)         2: (1, 1)
%   4: (2, -2)        3: (2, -2)
%   5: (2, -1)        
%   6: (2, 0)         4: (2, 0)     2: (2, 0)
%   7: ...
%   
%   For a demonstration on how to use this function,
%   see also BASISFUNCTIONSTEST, INDEXCONVERSIONTEST
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

switch spatialDimension
    case 1
        l = linIdx;
        m = zeros(size(l));
    case 2
        l = ceil(-3 / 2 + sqrt( 9 / 4 + 2 * linIdx) - 1e-15);
%     l = floor(-1 / 2 + sqrt( 1 / 4 + 2 * linIdx) + 1e-15);
        m = 2 * (linIdx - l .* (l + 1) / 2) - l;
    case 3
        l = floor(sqrt(linIdx) + 1e-15);
        m = linIdx - l.^2 - l;
    otherwise
        error('Invalid spatial dimension!')
end
end