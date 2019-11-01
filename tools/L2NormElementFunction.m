function E = L2NormElementFunction(fAtElements, points, connectivityList)
% L2NORMELEMENTFUNCTION Compute the L2 norm of a piecewise-constant
%   function w.r.t. the elements of a triangular mesh.
% 
%   For a demonstration on how to use this function,
%   see also TOOLSTEST
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%
if size(connectivityList, 1) == 3 && isvector(fAtElements)...
        && length(fAtElements) == size(connectivityList, 2)
    
    elementVol = getElementVolumes(points, connectivityList);
    E = sqrt(sum((fAtElements(:)').^2 .* elementVol));
else
   error('Dimension error in connectivity list!') 
end
end