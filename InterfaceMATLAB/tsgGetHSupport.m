function [coefficients] = tsgGetHSupport(lGrid)
%
% [coefficients] = tsgGetHSupport(lGrid)
%
% retrieves the supprot of the hierarchical basis functions
%
% INPUT:
%
% lGrid: a grid list created by tsgMake***(...) command
%
% OUTPUT:
%
% coefficients: the support of the basis functions in an array of size [num_points, num_dimensions]
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC] = tsgMakeFilenames(lGrid.sName);

sCommand = [sTasGrid,' -gethsupport'];

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
