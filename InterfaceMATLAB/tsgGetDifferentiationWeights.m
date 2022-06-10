function [weights] = tsgGetDifferentiationWeights(lGrid, mX)
%
% [weights] = tsgGetDifferentiationWeights(lGrid, mX)
%
% it gives the weights for differentiation weights
%
% INPUT:
%
% lGrid: a grid list created by tsgMakeXXX(...) command
%
% mX: an array of size [num_x, dimensions]
%     specifies the points where the derivative should be evaluated
%     Note: do not confuse points here with the nodes of the grid
%           here points are user specified points to evaluate the
%           derivative
%
% OUTPUT:
%
% weights: an array of size [num_x, number_of_points]
%          the values associated with the derivative at those points
%
%  The Jacobian of f(mX(i,:)) is approximated by:
%
%    f(points)' * weights(:, (i-1)*dimensions+1:i*dimensions)
%
%  where points is obtained from the tsgGetPoints() function.
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC] = tsgMakeFilenames(lGrid.sName);

sCommand = [sTasGrid,' -getdiffweights'];

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
