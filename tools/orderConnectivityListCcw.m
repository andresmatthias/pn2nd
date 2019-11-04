function [connectivityListOrdered, changeFlag] = orderConnectivityListCcw(points, connectivityList)
% ORDER_CONNECTIVITY_LIST_CCW Order the vertices of each triangle of a 
%   triangular mesh counter clockwise.
%
%   For a demonstration on how to use this function,
%   see also MESHAUXTEST
%
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

    changeFlag = 0;
    if size(connectivityList, 1) ~= 3
        error('Connectivity list not shaped properly!')
    end
    if size(points, 1) ~= 2
        error('Connectivity list not shaped properly!')
    end
    tol     = 1e-14;
    connectivityListOrdered  = connectivityList;
    for t = 1 : size(connectivityList, 2)
        A       = points(:, connectivityList(1, t));
        B       = points(:, connectivityList(2, t));
        C       = points(:, connectivityList(3, t));
        dAB     = B - A;
        dAB     = dAB / norm(dAB);
        dBC     = C - B;
        dBC     = dBC / norm(dBC);
        comp3   = dAB(1) * dBC(2) - dAB(2) * dBC(1);
        if abs(comp3) < tol
            warning('Vertices of element %1.0f close to collinear: cross_3 = %e', t, comp3)
        elseif  comp3 < 0
            tmp         = connectivityListOrdered(2, t);
            connectivityListOrdered(2, t) = connectivityListOrdered(3, t);
            connectivityListOrdered(3, t) = tmp;
            changeFlag = 1;
        end
    end
end