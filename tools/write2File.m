function write2File(fileName, headers, matrix, varargin)
    if nargin > 3
        prec = varargin{1};
    else
       prec = 5;
    end
    fid = fopen(fileName, 'wt');
    tmp = sprintf('%s,', headers{:});
    tmp = tmp(1 : end - 1);
    fprintf(fid, tmp);
    fprintf(fid, '\n');
    dlmwrite(fileName, matrix, '-append', 'precision', prec);
    fclose(fid);
end