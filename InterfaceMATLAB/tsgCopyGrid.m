function  [lNewGrid] = tsgCopyGrid(lOldGrid, sNewGridName)
%
% [lNewGrid] = tsgCopyGrid(lOldGrid, sNewGridName)
%
%  Makes a physical copy of a grid that differs only in the name
%  Note: tsgCopyGrid will create a new grid file in the work folder, while
%  the command lNewGrid = lOldGrid will only create an alias
%
% INPUT:
%
% lGrid: a grid list created by tsgMake***(...)
%
% sNewGridName: the name for the new grid, the new grid will be identical
% to the old one in every other way
%
% OUTPUT:
%
% lNewGrid: a grid object pointing to the physical copy of the old grid
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC] = tsgMakeFilenames(lOldGrid.sName);

lNewGrid = lOldGrid;
lNewGrid.sName = sNewGridName;
[sFileGNew, sFileX, sFileV, sFileO, sFileW, sFileC, sFileL] = tsgMakeFilenames(lNewGrid.sName);

sFileG    = regexprep(sFileG, '\\ ', ' ');
sFileGNew = regexprep(sFileGNew, '\\ ', ' ');
[status, cmdout] = copyfile(sFileG, sFileGNew);

if (~isempty(cmdout))
    fprintf(1,['WARNING: Command had non-empty output:\n']);
    disp(cmdout);
end

end
