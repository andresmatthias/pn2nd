function linIdx = degOrder2linearIdx(l, m, spatialDimension)
% DEGORDER2LINEARIDX Convert a tuple with degree and order of the 
%   corresponding real spherical harmonic (S_l^m: l:degree, m: order) to
%   the linear index as a basis function (when serially numbering all
%   basis functions, depending on symmetry assumptions (spatialDimension)).
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
%   see also INDEXCONVERSIONTEST
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

switch spatialDimension
    case 1
        linIdx = l;
    case 2
        linIdx = l .* (l + 1) / 2 + (m + l) / 2;
    case 3
        linIdx = l.^2 + m + l;
    otherwise
        error('Invalid spatial dimension!')
end    
end