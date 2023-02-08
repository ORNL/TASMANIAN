function tsgLoadConstructedPoints(lGrid, mX, mY)
%
%   tsgLoadConstructedPoints(lGrid, mX, mY)
%
% loads the combination of inputs (mX) and outputs (mY) into the grid
% the points don't need to be in any particular order, but
% the i-th row of mX has to correspond to the i-th row of mY
% the points must be aligned to the gird, but may come from
% any level, e.g., grid level 3 can accept points from level 4 or up
% see loadConstructedPoints() in the on-line manual
%
% NOTE: this will call the C++ method beginConstruction() which will
%       clear any currently set refinement
%
% INPUT:
%
% lGrid: a grid list created by tsgMakeXXX(...)
%
% mX : matrix of size N by number of dimensions
%
% mY : matrix of size N by number of outputs
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC, sFileL] = tsgMakeFilenames(lGrid);

sCommand = [sTasGrid,' -loadconstructed'];

sCommand = [sCommand, ' -gridfile ', sFileG];

tsgWriteMatrix(sFileX, mX);

sCommand = [sCommand, ' -xf ', sFileX];
lClean.sFileX = 1;

tsgWriteMatrix(sFileV, mY);

sCommand = [sCommand, ' -vf ', sFileV];
lClean.sFileV = 1;

[status, cmdout] = system(sCommand);

if (max(size(strfind(cmdout, 'ERROR'))) ~= 0)
    disp(cmdout);
    error('The tasgrid execurable returned an error, see above');
    return;
else
    if (~ isempty(cmdout))
        fprintf(1,['Warning: Command had non-empty output:\n']);
        disp(cmdout);
    end
end

if (exist('lClean'))
    tsgCleanTempFiles(lGrid, lClean);
end

end
