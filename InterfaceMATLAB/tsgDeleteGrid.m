function tsgDeleteGrid(lGrid)
%
% tsgDeleteGrid(lGrid)
%
% deletes all of the background files used by the grid
%
% INPUT:
%
% lGrid: a grid list created by tsgMakeGrid(...)
%
% NOTE: lGrid gets deleted and it can no longer be used
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC] = tsgMakeFilenames(lGrid.sName);

sFileG = regexprep(sFileG, '\\ ', ' ');
if (exist(sFileG, 'file') == 2)
    delete(sFileG);
end
lFiles.all = 1;
tsgCleanTempFiles(lGrid, lFiles);

end
