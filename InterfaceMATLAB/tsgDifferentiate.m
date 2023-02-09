function [result] = tsgDifferentiate(lGrid, mX)
%
% [result] = tsgDifferentiate(lGrid, points)
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
% result: if num_x is 1, this is the [outputs, dimensions] sized
%         Jacobian array/matrix of the surrogate model at mX;
%         if num_x is >1, this is a cell array of size [num_x] whose i-th
%         entry is the [outputs, dimensions] sized Jacobian array/matrix
%         of the surrogate model at mX(i,:).
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC] = tsgMakeFilenames(lGrid);

sCommand = [sTasGrid,' -differentiate'];

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
    if (~ isempty(cmdout))
        fprintf(1, ['Warning: Command had non-empty output:\n']);
        disp(cmdout);
    end
    raw_result = tsgReadMatrix(sFileO);
    num_x = size(mX, 1);
    num_dims = size(mX, 2);
    num_outs = size(raw_result, 2) / num_dims;
    % NOTE: C++ TASMANIAN arrays are stored row-major, but MATLAB assumes 1D vectors are column-major.
    if (num_x == 1)
        result = reshape(raw_result, [num_dims, num_outs])';
    else
        result = cell(num_x, 1);
        for i=1:num_x
            result{i} = reshape(raw_result(i,:), [num_dims, num_outs])';
        end
    end
end

tsgCleanTempFiles(lGrid, lClean);

end
