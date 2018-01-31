function [coefficients] = tsgGetHCoefficients(lGrid)
%
% [coefficients] = tsgGetHCoefficients(lGrid)
%
% retrieves the get the hierarchical coefficients of the grid
%
% INPUT:
%
% lGrid: a grid list created by tsgMake***(...) command
%
% OUTPUT:
%
% coefficients: the hierarchical coefficients of the grid in an array of dimension [num_poits, iOut]
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC] = tsgMakeFilenames(lGrid.sName);

sCommand = [sTasGrid,' -getcoefficients'];

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
    [coefficients] = tsgReadMatrix(sFileO);
end

tsgCleanTempFiles(lGrid, lClean);

end
