function [surpluses] = tsgGetSurpluses(lGrid)
%
% [surpluses] = tsgGetSurpluses(lGrid)
%
% retrieves the surpluses associated with an existing grid
%
% INPUT:
%
% lGrid: a grid list created by tsgMakeXXX(...) command
%
% OUTPUT:
%
% surpluses: the surplus associated with each point on the grid
%            in an array of dimension [num_poits, iOut]
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC] = tsgMakeFilenames(lGrid.sName);

sCommand = [sTasGrid,' -getsurpluses'];

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
    [surpluses] = tsgReadMatrix(sFileO);
end

tsgCleanTempFiles(lGrid, lClean);

end