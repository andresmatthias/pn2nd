function createFolder(buildFolder)
% CREATEFOLDER Create folder if it doesn't exist already.
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

    if ~exist(buildFolder, 'dir')
        mkdir(buildFolder)
    end
end