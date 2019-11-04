function E = L2NormElementFunction(fAtElements, points, connectivityList)
% L2NORMELEMENTFUNCTION Compute the L2 norm of a piecewise-constant
%   function w.r.t. the elements of a triangular mesh.
% 
%   For a demonstration on how to use this function,
%   see also TOOLSTEST
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%
if size(connectivityList, 1) == 3 && isvector(fAtElements)...
        && length(fAtElements) == size(connectivityList, 2)
    
    elementVol = getElementVolumes(points, connectivityList);
    E = sqrt(sum((fAtElements(:)').^2 .* elementVol));
else
   error('Dimension error in connectivity list!') 
end
end