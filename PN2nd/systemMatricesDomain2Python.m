function [] = systemMatricesDomain2Python(N_PN, spatialDimension, kernelName,...
    buildFolder, systemMatricesDomain)
% SYSTEMMATRICESDOMAIN2PYTHON Write system matrices for the PDE system of 
%   the second-order formulation on the domain to Python file.
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

printPrecision = 8;
tab = '    ';

%% system matrices

fileID = fopen(sprintf('%ssystemMatricesDomain_d_%d_k_%s_N_%d.py', buildFolder,...
        spatialDimension, kernelName, N_PN), 'w');
fprintf(fileID, 'import fenics as fe \nimport numpy as np\n\n\ndef getSystemMatricesDomain(par):\n');

%% build file
parameterRenaming

fprintf(fileID, '\n#-------- system matrices --------\n');
fprintf(fileID, '\n%ssystemMatricesDomain = dict()', tab);
fn = fieldnames(systemMatricesDomain);
for k = 1:length(fn)
    printMatrix(fn{k}, systemMatricesDomain.(fn{k}));
end

fprintf(fileID, '\n%sreturn systemMatricesDomain', tab);

fclose(fileID);

%%
function [] = parameterRenaming()
    fprintf(fileID, '\n#-------- parameter renaming --------\n');
    fprintf(fileID, '%ssigma_a = par[''sigma_a'']\n', tab);
    fprintf(fileID, '%ssigma_s = par[''sigma_s'']\n', tab);
    fprintf(fileID, '%ssigma_t = fe.Expression(''sigma_a + sigma_s'', degree=1, sigma_a=sigma_a, sigma_s=sigma_s)\n', tab);
end

function [] = printMatrix(name, M)
    fprintf(fileID, '\n');
    aux = sprintf('%ssystemMatricesDomain[''%s''] = [', tab, name);
    for i = 1:size(M, 1)
       tmp = char(vpa(M(i, 1), printPrecision));
       tmp = strrep(tmp, '^','**');
       aux = [aux, '[', tmp];
        for j = 2:size(M, 2)
            tmp = char(vpa(M(i, j), printPrecision));
            tmp = strrep(tmp, '^','**');
            aux = [aux, ', ', tmp];
        end
        if i < size(M, 1)
            fprintf(fileID, sprintf('%s],\n', aux));
            aux = [tab, tab];
        end
    end
    fprintf(fileID, '%s]]\n', aux);
end

end