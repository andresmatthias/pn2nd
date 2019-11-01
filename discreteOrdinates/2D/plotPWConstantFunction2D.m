function plotPWConstantFunction2D(MSH, QP_values, flat)  
% PLOTPWCONSTANTFUNCTION2D Plot a piecewise constant function on a
%   triangular mesh.
%
%   If the function u is given by its Nodal values, use
%   QP_values = u(mesh.connectivityList), or barycentric interpolation.
%   If the function u is piecewise constant w.r.t. the elements, then use
%   QP_values = repmat(u, 1, 3)
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%
    quadrature_points = [0, 0; 1, 0; 0, 1]';
    MSH.points = MSH.points';
    MSH.connectivityList = MSH.connectivityList';
    plotAux(MSH, quadrature_points, QP_values, flat)   
   
end

function plotAux(MSH, quadrature_points,QP_values,flat)
    if nargin<3
        flat = true;
    end
    nx = size(MSH.connectivityList, 1);
    %plots values, format n_cells x n_quadrature points
    Triangles = MSH.connectivityList;
    P = MSH.points;
    [Qp,ia] = unique(quadrature_points','rows');
    QP_values = QP_values(:,ia);
    Tr_h = delaunayTriangulation(Qp);
    n_qp = size(Qp,1);
    n_tr = size(Tr_h.ConnectivityList,1);
    n_cells = nx(1);
    Triangles = Triangles(1:n_cells,:);
    P = P(1:max(Triangles(:)),:);
    Points_g=zeros(numel(QP_values),2);
    %CL_g=zeros(n_tr*n_cells,3);
    TriList = Tr_h.ConnectivityList;
    %for j=1:size(QP_values,2)
    CL_g = reshape(bsxfun(@plus,permute(TriList,[1 3 2]), n_qp*(0:n_cells-1)), n_tr*n_cells, 3);
    for j = 1:n_cells
       A = P(Triangles(j,1),:);
       B =P (Triangles(j,2),:);
       C =P (Triangles(j,3),:);
       X = [B-A; C-A];
       QP_transf = bsxfun(@plus,Qp*X,A);
       Points_g((j-1)*n_qp+1:(j-1)*n_qp+n_qp,:) = QP_transf;
    %CL_g(n_tr*(j-1)+1:n_tr*(j-1)+n_tr,:)=TriList+n_qp*(j-1);
    end
    % trisurf(CL_g,Points_g(:,1),Points_g(:,2),[QP_values(:)]);
    a = squeeze(QP_values)';
    
    cmap = viridis(256);
    colormap(cmap);
    if flat
        trisurf(CL_g, Points_g(:,1), Points_g(:,2),0*a(:),a(:));
    else
       trisurf(CL_g, Points_g(:,1), Points_g(:,2),[a(:)]);
    end
end 