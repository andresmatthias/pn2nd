function E = L2DifferenceMeshFunctionVsPwConstant(fPwLinear, fPwConst, points, connectivityList)
% L2DIFFERENCEMESHFUNCTIONVSPWCONSTANT Compute the L2 difference between a
%   piecewise linear function w.r.t. the elements of the triangular mesh 
%   (fPwLinear; given by its values on each node) and a piecewise constant 
%   function w.r.t. the elements of the triangular mesh (fConst; given by
%   its values on each element.
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

if size(connectivityList, 1) == 3 ...
        && isvector(fPwLinear) && isvector(fPwConst) ...
        && length(fPwConst) == size(connectivityList, 2)...
        && length(fPwLinear) == size(points, 2)
        
    fPwLinear = fPwLinear(:)';
    fPwConst = fPwConst(:)';
    
    elementVol = getElementVolumes(points, connectivityList);
    fAtA = fPwLinear(connectivityList(1, :)) - fPwConst;
    fAtB = fPwLinear(connectivityList(2, :)) - fPwConst;
    fAtC = fPwLinear(connectivityList(3, :)) - fPwConst;
    
    fAtMidAB = 1/2 * (fAtA + fAtB);
    fAtMidBC = 1/2 * (fAtB + fAtC);
    fAtMidAC = 1/2 * (fAtA + fAtC);
    
    E = sqrt(sum((1/3 * fAtMidAB.^2 + 1/3 * fAtMidBC.^2 + 1/3 * fAtMidAC.^2) .* elementVol));
else
   error('Dimensions are not correct.') 
end
end

