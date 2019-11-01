function E = L2NormMeshFunction(fAtPoints, points, connectivityList)
% L2NORMMESHFUNCTION Compute the L2 norm of a piecewise-linear
%   function w.r.t. the elements of a triangular mesh, given by its nodal
%   values.
%   
%   Quadrature rule on unit triangle (A=(0,0), B=(1,0), C=(0,1)) for
%   polynomials up to degree 2 (see, e.g., Numerical Treatment of Partial 
%   Differential Equations; Grossmann, Roos, Stynes): 
%   I = 1/6 * g(1/2, 1/2) + 1/6 * g(0, 1/2) + 1/6 * g(1/2, 0)
% 
%   For a demonstration on how to use this function,
%   see also TOOLSTEST
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

if size(connectivityList, 1) == 3 && size(points, 1) == 2 ...
    && isvector(fAtPoints) && length(fAtPoints) == size(points, 2)

    elementVol = getElementVolumes(points, connectivityList);
    fAtA = fAtPoints(connectivityList(1, :));
    fAtB = fAtPoints(connectivityList(2, :));
    fAtC = fAtPoints(connectivityList(3, :));
    
    fAtMidAB = 1/2 * (fAtA + fAtB);
    fAtMidBC = 1/2 * (fAtB + fAtC);
    fAtMidAC = 1/2 * (fAtA + fAtC);
    
    E = sqrt(sum((1/3 * fAtMidAB.^2 + 1/3 * fAtMidBC.^2 + 1/3 * fAtMidAC.^2) .* elementVol));
else
   error('Dimension error in connectivity list!') 
end
end