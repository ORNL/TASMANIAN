function [new_points] = tsgRefineAnisotropic(lGrid, sType, iMinNew, iOut)
%
% [new_points] = tsgRefineAnisotropic(lGrid, sType, iMinNew, iOut)
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
% OUTPUT:
%
% new_points: (optional) array
%             the new set of needed points
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC] = tsgMakeFilenames(lGrid.sName);

sCommand = [sTasGrid,' -refineaniso'];

sCommand = [sCommand, ' -gridfile ', sFileG];

sCommand = [sCommand, ' -type ', sType];
sCommand = [sCommand, ' -ming ', num2str(iMinNew)];

if (exist('iOut'))
    sCommand = [sCommand, ' -refout ', num2str(iOut)];
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