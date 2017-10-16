function [weights] = tsgGetInterpolationWeights(lGrid, mX)
%
% [weights] = tsgGetInterpolationWeights(lGrid, mX)
%
% it gives the weights for interpolation (or approximation)
%
% INPUT:
%
% lGrid: a grid list created by tsgMakeXXX(...) command
%
% mX: an array of size [num_x, dimensions]
%     specifies the points where the interpolant should be evaluated
%     Note: do not confuse points here with the nodes of the grid
%               here points are user specified points to evaluate the
%               interpolant (or approximation)
%
% OUTPUT:
%
% weights: an array of size [num_x, number_of_points]
%          the values associated with the interpolant at those points
%
%  f(mX(i,:)) is approximated by weights(i,:)' * f(points)
%  where points is obtained from the tsgGetPoints() function
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC] = tsgMakeFilenames(lGrid.sName);

sCommand = [sTasGrid,' -getinterweights'];

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
    if (~isempty(cmdout))
        fprintf(1,['Warning: Command had non-empty output:\n']);
        disp(cmdout);
    end
    [weights] = tsgReadMatrix(sFileO);
end

tsgCleanTempFiles(lGrid, lClean);

end
