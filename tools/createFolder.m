function createFolder(buildFolder)
% CREATEFOLDER Create folder if it doesn't exist already.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

    if ~exist(buildFolder, 'dir')
        mkdir(buildFolder)
    end
end