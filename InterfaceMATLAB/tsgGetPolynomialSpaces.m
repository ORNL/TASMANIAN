function [poly] = tsgGetPolynomialSpaces(lGrid, bInterpolate)
%
% [poly] = tsgGetPolynomialSpaces(lGrid, bInterpolate)
%
% returns a matrix corresponding to the polynomial space that is integrated
% or interpolated exactly
%
% NOTE: this can be called only for global and sequence grids
%
% INPUT:
%
% lGrid: a grid list created by tsgMakeXXX(...)
%
% bInterpolate: (boolean)
%               specifies whether to consider integration or interpolation
%
% OUTPUT:
%
% poly: (matrix of size num_polynomial_basis_functions X iDim)
%       The polynomial space is
%       span{  x.^poly(i,:)  }, for x in R^d and i = 1 ... size(poly, 1)
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC] = tsgMakeFilenames(lGrid.sName);

sCommand = [sTasGrid,' -getpoly'];

sCommand = [sCommand, ' -gridfile ', sFileG];

if (bInterpolate)
    sCommand = [sCommand, ' -type iptotal'];
else
    sCommand = [sCommand, ' -type qptotal'];
end

sCommand = [sCommand, ' -of ',sFileO];
lClean.sFileO = 1;

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
    [poly] = tsgReadMatrix(sFileO);
end

tsgCleanTempFiles(lGrid, lClean);

end