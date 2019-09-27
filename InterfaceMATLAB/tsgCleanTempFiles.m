function tsgCleanTempFiles(lGrid, lFiles)
%
% tsgCleanTempFiles(lGrid, bFlags)
%
% Cleans the temporary files associated with lGrid, i.e., all files other
% than the _FileG. If lFlags is [], then all files are cleaned. Otherwise
% only the files with entries in lFlags are cleaned, this is done to limit
% unecessary calls to the 'rm' command.
%
% You may call this function, but it exists mostly to be called by other
% functions. In fact, other functions should clean their temp files so you
% should never have to worry about this.
%
% INPUT:
%       lGrid:  a grid list created by a tsgMake***(...) command (i.e., an
%               existing grid)
%
%       lFiles: object with fields sFileX, sFileV, sFileO, sFileW, sFileC
%               corresponding to the files to be deleted. If the field
%               exist, then the file is deleted. If lFiles is omitted, then
%               all files are deleted (except the permanent grid file FileG
%

% generate filenames
[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC, sFileL] = tsgMakeFilenames(lGrid.sName);

sFileX = regexprep(sFileX, '\\ ', ' ');
sFileV = regexprep(sFileV, '\\ ', ' ');
sFileO = regexprep(sFileO, '\\ ', ' ');
sFileW = regexprep(sFileW, '\\ ', ' ');
sFileC = regexprep(sFileC, '\\ ', ' ');
sFileL = regexprep(sFileL, '\\ ', ' ');

if (isfield(lFiles, 'all'))
    if (exist(sFileX, 'file') == 2)
        delete(sFileX);
    end
    if (exist(sFileV, 'file') == 2)
        delete(sFileV);
    end
    if (exist(sFileO, 'file') == 2)
        delete(sFileO);
    end
    if (exist(sFileW, 'file') == 2)
        delete(sFileW);
    end
    if (exist(sFileC, 'file') == 2)
        delete(sFileC);
    end
    if (exist(sFileL, 'file') == 2)
        delete(sFileL);
    end
else
    if (isfield(lFiles, 'sFileX'))
        if (exist(sFileX, 'file') == 2)
            delete(sFileX);
        end
    end
    if (isfield(lFiles, 'sFileV'))
        if (exist(sFileV, 'file') == 2)
            delete(sFileV);
        end
    end
    if (isfield(lFiles, 'sFileO'))
        if (exist(sFileO, 'file') == 2)
            delete(sFileO);
        end
    end
    if (isfield(lFiles, 'sFileW'))
        if (exist(sFileW, 'file') == 2)
            delete(sFileW);
        end
    end
    if (isfield(lFiles, 'sFileC'))
        if (exist(sFileC, 'file') == 2)
            delete(sFileC);
        end
    end
    if (isfield(lFiles, 'sFileL'))
        if (exist(sFileL, 'file') == 2)
            delete(sFileL);
        end
    end
end

end
