function [new_points] = tsgRefineAnisotropic(lGrid, sType, iMinNew, iOut, vLimitLevels)
%
% [new_points] = tsgRefineAnisotropic(lGrid, sType, iMinNew,
%                                     iOut, vLimitLevels)
%
% computes anisotropic weights corresponding to sType and refines the grid
% so that at least iMinNew points are added to the grid while respecting
% the computed anisotropic weights
%
% NOTE: this can be called only for global and sequence grids
%
% INPUT:
%
% lGrid: a grid list created by tsgMakeXXX(...)
%
% sType: (string giving the tensor selection strategy)
%        'level'       'curved'         'tensor'         'iptensor'
%        'iptotal'     'ipcurved'       'qptotal'        'qpcurved'
%
% iMinNew: (integer strictly positive)
%          minimum number of new points to include in the new grid
%
% iOut: (integer giving the output to be used for the refinement)
%       selects which output to use for refinement
%
% vLimitLevels: (optional vector of integers of size iDim)
%               limit the level in each direction, no points beyond the
%               specified limit will be used, e.g., in 2D using
%               clenshaw-curtis rule, [1, 99] forces the grid to have
%               at most 3 possible values in the first variable and
%               ~2^99 (practicallyt infinite) number in the second
%               direction. vLimitLevels overwrites iMinNew, if using
%               vLimitLevels the number of new points may be less than
%               iMinNew, increase iMinNew until a desired number of
%               points is selected
%
% OUTPUT:
%
% new_points: (optional) array
%             the new set of needed points
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC, sFileL] = tsgMakeFilenames(lGrid.sName);

sCommand = [sTasGrid,' -refineaniso'];

sCommand = [sCommand, ' -gridfile ', sFileG];

sCommand = [sCommand, ' -type ', sType];
sCommand = [sCommand, ' -ming ', num2str(iMinNew)];

if (exist('iOut'))
    sCommand = [sCommand, ' -refout ', num2str(iOut)];
end

% set level limits
if (exist('vLimitLevels') && (max(size(vLimitLevels)) ~= 0))
    if (min(size(vLimitLevels)) ~= 1)
        error(' vLimitLevels must be a vector, i.e., one row or one column');
    end
    if (max(size(vLimitLevels)) ~= lGrid.iDim)
        error(' vLimitLevels must be a vector of size iDim');
    end
    if (size(vLimitLevels, 1) > size(vLimitLevels, 2))
        tsgWriteMatrix(sFileL, vLimitLevels');
    else
        tsgWriteMatrix(sFileL, vLimitLevels);
    end
    lClean.sFileL = 1;
    sCommand = [sCommand, ' -levellimitsfile ', sFileL];
end

% read the points for the grid
if (nargout == 1)
    sCommand = [sCommand, ' -of ',sFileO];
    lClean.sFileO = 1;
end

[status, cmdout] = system(sCommand);

if (max(size(findstr('ERROR', cmdout))) ~= 0)
    disp(cmdout);
    error('The tasgrid execurable returned an error, see above');
    return;
else
    if (~isempty(cmdout))
        fprintf(1,['WARNING: Command had non-empty output:\n']);
        disp(cmdout);
    end
    if (nargout == 1)
        [new_points] = tsgReadMatrix(sFileO);
    end
end

if (exist('lClean'))
    tsgCleanTempFiles(lGrid, lClean);
end

end
