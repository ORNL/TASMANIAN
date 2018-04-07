function [result] = tsgEvaluate(lGrid, mX)
%
% [result] = tsgEvaluate(lGrid, points)
%
% evaluates the intepolant at the points of interest and returns the result
% this should be called after the grid has been created and after values
% have been loaded
%
% INPUT:
%
% lGrid: a grid list created by tsgMakeXXX(...)
%
% mX: an array of size [num_x, dimensions]
%     specifies the points where the interpolant should be evaluated
%     Note: do not confuse points here with the nodes of the grid
%               here points are user specified points to evaluate the
%               interpolant (or approximation)
%
% NOTE: if lGrid has a field gpuDevice and if Tasmanian is build with
%       CUBLAS or CUDA options, then this will attempt to use the GPU
%       for acceleration. The gpuDevice should be an integer corresponding
%       to a valid Nvidia CUDA device, run tsgCoreTests() to see the
%       device list visible to Tasmanian and the corresponding number
%
% OUTPUT:
%
% result: an array of size [num_x, iOut]
%         the values of the interpolant at the corresponding points
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC] = tsgMakeFilenames(lGrid.sName);

sCommand = [sTasGrid,' -evaluate'];

sCommand = [sCommand, ' -gridfile ', sFileG];

tsgWriteMatrix(sFileX, mX);

sCommand = [sCommand, ' -xf ', sFileX];
lClean.sFileX = 1;

sCommand = [sCommand, ' -of ', sFileO];
lClean.sFileO = 1;

if (isfield(lGrid, 'gpuDevice'))
    sCommand = [sCommand, ' -gpuid ', num2str(lGrid.gpuDevice)];
end

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
