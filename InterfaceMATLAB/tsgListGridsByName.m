function tsgListGridsByName()
%
% tsgListGridsByName()
%
% prints a list of all grids in the WorkFiles folder regardless if they are
% associated with lGrid class
%

[sFiles, sTasGrid] = tsgGetPaths();

lFileList = dir(sFiles);

lGridIndexes = [];

for iI = 1:length(lFileList)
    if (~isempty(findstr('_FileG', lFileList(iI).name)))
        lGridIndexes = [lGridIndexes, iI];
    end
end

disp([' Grids in the work folder:']);
disp([' ']);

for iI = 1:length(lGridIndexes)
    line = lFileList(lGridIndexes(iI)).name;
    disp(line(1:end-6));
end

disp([' ']);

end