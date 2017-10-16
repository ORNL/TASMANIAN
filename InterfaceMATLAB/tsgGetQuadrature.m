function [weights, points] = tsgGetQuadrature(lGrid)
%
% [weights, points] = tsgGetQuadrature(lGrid)
%
% retrieves the quadrature points and weights from an existing grid
%
% INPUT:
%
% lGrid: a grid list created by a tsgMakeXXX(...) command
%
% OUTPUT:
%
% weights: the quadrature weights
%
% points: the quadrature nodes, i.e., quad(f) = weights' * f(points)
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC] = tsgMakeFilenames(lGrid.sName);

sCommand = [sTasGrid,' -getquadrature'];

sCommand = [sCommand, ' -gridfile ', sFileG];

sCommand = [sCommand, ' -of ', sFileO];
lClean.sFileO = 1;

[status, cmdout] = system(sCommand);

if (size(findstr('ERROR', cmdout)) ~= [0, 0])
    disp(cmdout);
    error('The tasgrid execurable returned an error, see above');
    return;
else
    if (~ isempty(cmdout))
        fprintf(1,['Warning: Command had non-empty output:\n']);
        disp(cmdout);
    end
    [wp] = tsgReadMatrix(sFileO);

    weights = wp(:,1);
    points = wp(:,2:end);
end

tsgCleanTempFiles(lGrid, lClean);

end