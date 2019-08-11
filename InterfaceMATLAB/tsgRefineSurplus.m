function [new_points] = tsgRefineSurplus(lGrid, fTolerance, sRefinementType, iOut, vLimitLevels, mCorrection)
%
% [new_points] = tsgRefineSurplus(lGrid, fTolerance, sRefinementType,
%                                 iOut, vLimitLevels)
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
%                  'classic'  'parents'   'direction'   'fds' 'stable'
%                  only for Local Polynomial and Wavelet grids
%
% iOut: (integer giving the output to be used for the refinement)
%       selects which output to use for refinement, must be specified for
%       Global Grids, other grid can omit this or set it to [] to indicate
%       the use of all outputs
%
% vLimitLevels: (optional vector of integers of size iDim)
%               limit the level in each direction, no points beyond the
%               specified limit will be used, e.g., in 2D using
%               localp rule, [1, 99] forces the grid to have
%               at most 3 possible values in the first variable and
%               ~2^99 (practicallyt infinite) number in the second
%               direction.
%
% mCorrection: (optional matrix with size num_points by num_outputs or 1)
%              the surpluses in the refinement procedure will be multiplied
%              by the correction
%              if iOut is -1, then one weight must be given per output and
%              the columns must match the number of outputs
%              if iOut is specified, mCorrection must have one column
%
% OUTPUT:
%
% new_points: (optional) array
%             the new set of needed points
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC, sFileL] = tsgMakeFilenames(lGrid.sName);

sCommand = [sTasGrid,' -refinesurp'];

sCommand = [sCommand, ' -gridfile ', sFileG];

sCommand = [sCommand, ' -tolerance ', num2str(fTolerance,17)];

if ((exist('sRefinementType')) && (max(size(sRefinementType)) ~= 0))
    sCommand = [sCommand, ' -reftype ', sRefinementType];
end

if (exist('iOut') && (max(size(iOut)) ~= 0))
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

% set correction
if (exist('mCorrection') && (max(size(mCorrection)) ~= 0))
    if (min(size(mCorrection)) ~= 1)
        error(' mCorrection must be a vector, i.e., one row or one column');
    end
    tsgWriteMatrix(sFileV, mCorrection);
    lClean.sFileV = 1;
    sCommand = [sCommand, ' -valsfile ', sFileV];
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
