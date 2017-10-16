function [new_points] = tsgRefineSurplus(lGrid, fTolerance, sRefinementType, iOut)
%
% [new_points] = tsgRefineSurplus(lGrid, fTolerance, sRefinementType, iOut)
%
% adds new points to the grid in the neighbourhood of existing points
% associated with large hierarchical surplus
%
% This function works for Local Polynomial and Wavelet grids only 
%
% INPUT:
%
% lGrid: a grid list created by tsgMakeXXX(...)
%
% fTolerance: (real non-negative number)
%             refine only in neighborhood of points that satisfy
%             max(point_surplus / max_value) > fTolerance
%             where max(...) is taken over all outputs
%             point_surplus is the hierarchical surplus
%             max_value is the largest loaded value associated with each output
%
%
% sRefinementType: (string indicating the refinement type)
%                  'classic'  'parents'   'direction'   'fds'
%                  only for Local Polynomial and Wavelet grids
%
% iOut: (integer giving the output to be used for the refinement)
%       selects which output to use for refinement, only for Global lGrid
%
% OUTPUT:
%
% new_points: (optional) array
%             the new set of needed points
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC] = tsgMakeFilenames(lGrid.sName);

sCommand = [sTasGrid,' -refinesurp'];

sCommand = [sCommand, ' -gridfile ', sFileG];

sCommand = [sCommand, ' -tolerance ', num2str(fTolerance,17)];

if ((exist('sRefinementType')) && (max(size(sRefinementType)) ~= 0))
    sCommand = [sCommand, ' -reftype ', sRefinementType];
end

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