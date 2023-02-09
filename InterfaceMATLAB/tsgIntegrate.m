function [result] = tsgIntegrate(lGrid)
%
% [result] = tsgEvaluate(lGrid, points)
%
% returns the integral of the interpolant
%
% INPUT:
%
% lGrid: a grid list created by tsgMakeXXX(...)
%
% OUTPUT:
%
% result: a vector of size [iOut]
%         the values of the integrals of each output
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC] = tsgMakeFilenames(lGrid);

sCommand = [sTasGrid,' -integrate'];

sCommand = [sCommand, ' -gridfile ', sFileG];

sCommand = [sCommand, ' -of ', sFileO];
lClean.sFileO = 1;

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
    [result] = tsgReadMatrix(sFileO);
end

tsgCleanTempFiles(lGrid, lClean);

end
