function [points, connectivityList] = importGmsh(pathToMesh, varargin)
% IMPORTGMSH Extract triangulation (points, connectivity list) from Gmsh's
%   .msh format.
% 
%   For a demonstration on how to use this function,
%   see also MESHAUXTEST
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%
    fid = fopen(pathToMesh);

    tline = fgets(fid);
    p_string = cell(0,4);
    TR_string = cell(0,8);
    while ischar(tline)
        if strcmp(tline(1:end-1),'$Nodes') == 1
            tline = fgets(fid);% jump to next row
            tline = fgets(fid); % kill number of nodes
            while strcmp(tline(1:end-1),'$EndNodes') ==0
                p_string(end+1,:) = (strsplit(tline,' '));%tline(2:end);
                tline = fgets(fid);
            end
        end        
        tline = fgets(fid);
        if strcmp(tline(1:end-1),'$Elements') == 1
            tline = fgets(fid);% jump to next row
            tline = fgets(fid); % kill number of elements
            while strcmp(tline(1:end-1),'$EndElements') ==0
                tmp = strsplit(tline,' ');
                if strcmp(tmp(2),'2')
                    TR_string(end+1,:) = strsplit(tline,' '); %tline(4:end);
                    tline = fgets(fid);
                else
                    tline = fgets(fid);
                end
            end
        end
    end
    fclose(fid);
    
    points = cellfun(@str2num, p_string);
    connectivityList = cellfun(@str2num, TR_string);
    points = points(:,2:end)';
    connectivityList = connectivityList(:,end-2:end)';
    if strcmp(varargin{1}, '2D')
       points = points(1:2, :);
    end
end