function [result] = tsgDifferentiate(lGrid, mX)
%
% [result] = tsgEvaluate(lGrid, points)
%
% differentiates the interpolant at the points of interest and returns the result
% this should be called after the grid has been created and after values
% have been loaded
%
% INPUT:
%
% lGrid: a grid list created by tsgMakeXXX(...)
%
% mX: an array of size [num_x, dimensions]
%     specifies the points where the interpolant should be differentiated
%     Note: do not confuse points here with the nodes of the grid
%               here points are user specified points to differentiate the
%               interpolant (or approximation)
%
% OUTPUT:
%
% result: an array of size [num_x, iOut * dimensions]
%         the derivative of the interpolant at the corresponding points
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC] = tsgMakeFilenames(lGrid.sName);

sCommand = [sTasGrid,' -differentiate'];

sCommand = [sCommand, ' -gridfile ', sFileG];

tsgWriteMatrix(sFileX, mX);

sCommand = [sCommand, ' -xf ', sFileX];
lClean.sFileX = 1;

sCommand = [sCommand, ' -of ', sFileO];
lClean.sFileO = 1;

[status, cmdout] = system(sCommand);

if (size(findstr('ERROR', cmdout)) ~= [0, 0])
    disp(cmdout);
    error('The tasgrid execurable returned an error, see above');
    return;
else
    if (~ isempty(cmdout))
        fprintf(1, ['Warning: Command had non-empty output:\n']);
        disp(cmdout);
    end
    [result] = tsgReadMatrix(sFileO);
end

tsgCleanTempFiles(lGrid, lClean);

end
