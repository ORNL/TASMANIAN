function [points] = tsgGetPoints(lGrid)
%
% [points] = tsgGetPoints(lGrid)
%
% retrieves the points associated with an existing grid
%
% INPUT:
%
% lGrid: a grid list created by tsgMakeXXX(...) command
%
% OUTPUT:
%
% points: the points of the grid in an array of dimension [num_poits, iDim]
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC] = tsgMakeFilenames(lGrid.sName);

sCommand = [sTasGrid,' -getpoints'];

sCommand = [sCommand, ' -gridfile ', sFileG];

sCommand = [sCommand, ' -of ', sFileO];
lClean.sFileO = 1;

[status, cmdout] = system(sCommand);

if (size(findstr('ERROR', cmdout)) ~= [0, 0])
    disp(cmdout);
    error('The tasgrid execurable returned an error, see above');
    return;
else
    if (~ isempty(cmdout))
        fprintf(1,['Warning: Command had non-empty output:\n']);
        disp(cmdout);
    end
    [points] = tsgReadMatrix(sFileO);
end

tsgCleanTempFiles(lGrid, lClean);

end