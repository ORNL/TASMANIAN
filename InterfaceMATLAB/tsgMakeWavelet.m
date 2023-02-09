function [lGrid, points] = tsgMakeWavelet(sGridName, iDim, iOut, iDepth, iOrder, mTransformAB, sConformalMap, vConfromalWeights, vLimitLevels)
%
% [lGrid, points] = tsgMakeWavelet(sGridName, iDim, iOut, iDepth, iOrder,
%                        mTransformAB, sConformalMap, vConfromalWeights)
%
% creates a new sparse grid using a sequence rule
%
% INPUT:
%
% sGridName: the name of the grid, give it a string name,
%            i.e. 'myGrid' or '1' or 'pi314'
%            DO NOT LEAVE THIS EMPTY
%
% iDim: (integer, positive)
%       the number of inputs
%
% iOut: (integer, non-negative)
%       the number of outputs
%
% iDepth: (integer non-negative)
%          controls the density of the grid, i.e., the number of
%          levels to use
%
% iOrder: (integer must be 1 or 3)
%         note that only wavelets of order 1 and 3 are implemented
%
% mTransformAB: (optional matrix of size iDim x 2)
%               for all but gauss-laguerre and gauss-hermite grids, the
%               transform specifies the lower and upper bound of the domain
%               in each direction. For gauss-laguerre and gauss-hermite
%               grids, the transform gives the a and b parameters that
%               change the weight to
%               exp(-b (x - a))  and  exp(-b (x - a)^2)
%
% sConformalMap: (optional string giving the type of transform)
%                conformal maps provide a non-linear domain transform,
%                approximation (quadrature or interpolation) is done
%                on the composition of f and the transform. A suitable
%                transform could reduce the error by as much as an
%                order of magnitude.
%
%                'asin': truncated MacLaurin series of arch-sin
%
% vConfromalWeights: (optional parameters for the conformal trnasform)
%               'asin': indicate the number of terms to keep after
%                       truncation
%
% vLimitLevels: (optional vector of integers of size iDim)
%               limit the level in each direction, no points beyond the
%               specified limit will be used, e.g., in 2D [1, 99] forces
%               the grid to have at most 3 possible values in the first
%               variable and ~2^99 (practicallyt infinite) number in the
%               second direction. vLimitLevels works in conjunction with
%               iDepth, for each direction, we chose the lesser of the
%               vLimitLevels and iDepth
%
% OUTPUT:
%
% lGrid: list containing information about the sparse grid, can be used
%        to call other functions
%
% points: (optional) the points of the grid in an array
%         of dimension [num_poits, dim]
%
% [lGrid, points] = tsgMakeWavelet(sGridName, iDim, iOut, iDepth, iOrder,
%                        mTransformAB, sConformalMap, vConfromalWeights)
%

% create lGrid object
lGrid.sName = sGridName;
lGrid.sFilename = tsgMakeGridFilename(sGridName);
lGrid.iDim  = iDim;
lGrid.iOut  =  iOut;
lGrid.sType = 'wavelet';

% check for conflict with tsgMakeQuadrature
if (strcmp(sGridName, ''))
    error('sGridName cannot be empty');
end

% generate filenames
[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC, sFileL] = tsgMakeFilenames(lGrid);

sCommand = [sTasGrid,' -makewavelet'];

sCommand = [sCommand, ' -gridfile ',   sFileG];
sCommand = [sCommand, ' -dimensions ', num2str(lGrid.iDim)];
sCommand = [sCommand, ' -outputs ',    num2str(lGrid.iOut)];
sCommand = [sCommand, ' -depth ',      num2str(iDepth)];
sCommand = [sCommand, ' -order ',      num2str(iOrder)];

% set the domain transformation
if (exist('mTransformAB') && (max(size(mTransformAB)) ~= 0))
    if (size(mTransformAB, 2) ~= 2)
        error(' mTransformAB must be a matrix with 2 columns');
    end
    if (size(mTransformAB, 1) ~= lGrid.iDim)
        error(' mTransformAB must be a matrix with iDim number of rows');
    end
    tsgWriteMatrix(sFileV, mTransformAB);
    lGlean.sFileV = 1;
    sCommand = [sCommand, ' -tf ',sFileV];
end

% set conformal mapping
if (exist('sConformalMap')  && (max(size(sConformalMap)) ~= 0))
    if (~exist('vConfromalWeights'))
        error(' sConformalMap requires vConfromalWeights')
    end
    sCommand = [sCommand, ' -conformaltype ', sConformalMap];
    if (size(vConfromalWeights, 1) > size(vConfromalWeights, 2))
        tsgWriteMatrix(sFileC, vConfromalWeights');
    else
        tsgWriteMatrix(sFileC, vConfromalWeights);
    end
    lClean.sFileC = 1;
    sCommand = [sCommand, ' -conformalfile ',sFileC];
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
if (nargout > 1)
    sCommand = [sCommand, ' -of ',sFileO];
    lGlean.sFileO = 1;
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
    if (nargout > 1)
        points = tsgReadMatrix(sFileO);
    end
end

if (exist('lClean'))
    tsgCleanTempFiles(lGrid, lClean);
end

end
