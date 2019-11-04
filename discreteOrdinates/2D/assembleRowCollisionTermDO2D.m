function indexSetCollisionLhs = assembleRowCollisionTermDO2D(tupidx, par, mesh)
% ASSEMBLEROWCOLLISIONTERMDO2D Assemble rows of the discrete ordinates
%   system matrix regarding the collision term in the original kinetic
%   equation.
% 
% tupidx:  first row: element indices: 
%          second row: ordinate indices: 
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

%% handles the following part
%
% $$\frac{1}{|T_i|}\int_{T_i} -\sigma_t I\textup{d}x + \frac{1}{|T_i|}\int_{T_i} \int\kappa(\Omega,\Omega') I \textup{d}\Omega'\textup{d}x$$
%
nIdx = size(tupidx, 2);
nOrd = par.nDiscreteOrdinates;
elementNodeIDs = mesh.connectivityList(:, tupidx(1, :)); % indices of element's vertices
pointsAux = mesh.points(:, elementNodeIDs(:));
matrixEntries = zeros(nOrd, nIdx);

% this one is huge: nr_ordinates^2*nr_elmements
adjacencyThroughKernel = [reshape(repmat(tupidx(1, :), nOrd, 1), 1, []);
                            repmat(1 : nOrd, 1, nIdx)];
adjacentIndices = reshape(tup2linDO(adjacencyThroughKernel, nOrd), nOrd, nIdx);

%% first part
%
% $$\frac{1}{|T_i|}\int_{T_i} -\sigma_t I\textup{d}x $$
%

sigma_t_ElementNodes = reshape(par.sigma_t(pointsAux), 3, []);
sigma_t_mean = sum(sigma_t_ElementNodes, 1) / 3;
idx_diagonal = sub2ind(size(matrixEntries), tupidx(2, :), 1 : nIdx);
matrixEntries(idx_diagonal) = -sigma_t_mean;

%% second part
%
% $$ \frac{1}{|T_i|}\int_{T_i} \int\kappa(\Omega,\Omega') I \textup{d}\Omega'\textup{d}x$$
%
sigma_s_ElementNodes = reshape(par.sigma_s(pointsAux), 3, []);
sigma_s_mean = sum(sigma_s_ElementNodes, 1) / 3;
matrixEntries = matrixEntries + repmat(sigma_s_mean, nOrd, 1) .* ...
    (repmat(par.ordinates.weights', 1, nIdx) .* par.kernelMatrix(:, tupidx(2, :)));

%% create index set for sparse structure
% attention to minus!
indexSetCollisionLhs = [reshape(repmat(1 : nIdx, par.nDiscreteOrdinates, 1), [], 1),...
    adjacentIndices(:), -matrixEntries(:)];
end