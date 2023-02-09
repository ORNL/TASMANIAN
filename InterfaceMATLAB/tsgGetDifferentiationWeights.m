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
% weights: if num_x == 1, this is the [num_points, dimensions] sized array/matrix
%          of derivative weights for the point mX;
%          if num_x > 2, this is a cell array of size [num_x] whose i-th
%          entry is a [num_points, dimensions] sized array/vector of derivative
%          weights for mX(i,:);
%
%  The Jacobian of f(mX(i,:)) is approximated by:
%
%    f(points)' * weights,    if num_x == 1
%    f(points)' * weights{i}, if num_x >= 2
%
%  where points is obtained from the tsgGetPoints() function.
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC] = tsgMakeFilenames(lGrid);

sCommand = [sTasGrid,' -getdiffweights'];

sCommand = [sCommand, ' -gridfile ', sFileG];

tsgWriteMatrix(sFileX, mX);

sCommand = [sCommand, ' -xf ', sFileX];
lClean.sFileX = 1;

sCommand = [sCommand, ' -of ', sFileO];
lClean.sFileO = 1;

[status, cmdout] = system(sCommand);

if (max(size(strfind(cmdout, 'ERROR'))) ~= 0)
    disp(cmdout);
    error('The tasgrid execurable returned an error, see above');
    return;
else
    if (~isempty(cmdout))
        fprintf(1,['Warning: Command had non-empty output:\n']);
        disp(cmdout);
    end
    raw_result = tsgReadMatrix(sFileO);
    num_x = size(mX, 1);
    num_dims = size(mX, 2);
    num_points = size(raw_result, 2) / num_dims;
    % NOTE: C++ TASMANIAN arrays are stored row-major, but MATLAB assumes 1D vectors are column-major.
    if (num_x == 1)
        weights = reshape(raw_result, [num_dims, num_points])';
    else
        weights = cell(num_x, 1);
        for i=1:num_x
            weights{i} = reshape(raw_result(i,:), [num_dims, num_points])';
        end
    end
end

tsgCleanTempFiles(lGrid, lClean);

end
