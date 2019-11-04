function indexSetCollisionLhs = assembleRowCollisionTermDO1D(tupidx, par)
% ASSEMBLEROWCOLLISIONTERMDO1D Assemble rows of the discrete ordinates
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
nOrdinates = par.nDiscreteOrdinates;

matrixEntries = zeros(nOrdinates, nIdx);

% this one is huge: nr_ordinates^2*nr_elmements
adjacencyThroughKernel = [reshape(repmat(tupidx(1, :), nOrdinates, 1), 1, []);
                            repmat(1 : nOrdinates, 1, nIdx)];
adjacentIndices = reshape(tup2linDO(adjacencyThroughKernel, nOrdinates), nOrdinates, nIdx);

%% first part
%
% $$\frac{1}{|T_i|}\int_{T_i} -\sigma_t I\textup{d}x $$
%
sigma_t_mean = (par.sigma_t(par.grid(tupidx(1, :) + 1)) + par.sigma_t(par.grid(tupidx(1, :)))) / 2;

idxDiagonal = sub2ind(size(matrixEntries), tupidx(2, :), 1 : nIdx);
matrixEntries(idxDiagonal) = -sigma_t_mean;

%% second part
%
% $$ \frac{1}{|T_i|}\int_{T_i} \int\kappa(\Omega,\Omega') I \textup{d}\Omega'\textup{d}x$$
%
sigma_s_mean = (par.sigma_s(par.grid(tupidx(1, :) + 1)) + par.sigma_s(par.grid(tupidx(1, :)))) / 2;
matrixEntries = matrixEntries + repmat(sigma_s_mean, nOrdinates, 1) .* ...
    (repmat(par.ordinates.weights', 1, nIdx) .* par.kernelMatrix(:, tupidx(2, :)));

%% create index set for sparse structure
% attention to minus!
indexSetCollisionLhs = [reshape(repmat(1 : nIdx, par.nDiscreteOrdinates, 1), [], 1),...
    adjacentIndices(:), -matrixEntries(:)];
end
