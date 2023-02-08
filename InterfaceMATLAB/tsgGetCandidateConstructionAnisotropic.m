function [new_points] = tsgGetCandidateConstructionAnisotropic(lGrid, sType, iOutOrAnisotropy, vLimitLevels)
%
% [new_points] = tsgGetCandidateConstructionAnisotropic(lGrid,
%                              sType, iOutputOrAnisotropy, vLimitLevels)
%
% returns points for asynchronous (dynamic) construction where
% the "more important" points are sorted first, the sorting is done
% using anisotropic criteria
% see getCandidateConstructionPoints() in the on-line manual
%
% NOTE: this will call the C++ method beginConstruction() which will
%       clear any currently set refinement
%
% INPUT:
%
% lGrid: a grid list created by tsgMakeXXX(...)
%
% sType: (string giving the tensor selection strategy)
%        'level'       'curved'         'tensor'         'iptensor'
%        'iptotal'     'ipcurved'       'qptotal'        'qpcurved'
%
% iOutOrAnisotropy: (integer or vector of integers)
%          if it is a single integer, it will be used to indicate
%          the output to use to compute the ansiotropic weights
%          (could be -1 for Sequence and Fourier grids)
%          could be a vector of integers indicating anisotropic weights
%          which will skip computing the weights
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
% new_points: array
%             the new set of sorted needed points
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC, sFileL] = tsgMakeFilenames(lGrid);

sCommand = [sTasGrid,' -getconstructpnts'];

sCommand = [sCommand, ' -gridfile ', sFileG];

sCommand = [sCommand, ' -type ', sType];

if (max(size(iOutOrAnisotropy)) == 1)
    sCommand = [sCommand, ' -refout ', num2str(iOutOrAnisotropy)];
else
    if (min(size(iOutOrAnisotropy)) ~= 1)
        error(' iOutOrAnisotropy must be a vector, i.e., one row or one column');
    end
    if (max(size(iOutOrAnisotropy)) ~= lGrid.iDim)
        error(' iOutOrAnisotropy must be a vector of size iDim');
    end
    if (size(iOutOrAnisotropy, 1) > size(iOutOrAnisotropy, 2))
        tsgWriteMatrix(sFileW, iOutOrAnisotropy');
    else
        tsgWriteMatrix(sFileW, iOutOrAnisotropy);
    end
    lClean.sFileW = 1;
    sCommand = [sCommand, ' -anisotropyfile ', sFileW];
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

if (max(size(strfind(cmdout, 'ERROR'))) ~= 0)
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
