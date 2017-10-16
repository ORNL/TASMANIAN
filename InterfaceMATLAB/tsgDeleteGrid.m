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

if (exist(sFileG, 'file') == 2)
    delete(sFileG);
end
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

end