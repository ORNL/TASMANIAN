function tsgLoadHierarchical(lGrid, mValues)
%
% tsgLoadValues(lGrid, mValues)
%
% loads the values of the target function at the needed points
%
% INPUT:
%
% lGrid: a grid list created by tsgMakeXXX(...) command
%
% mValues: a matrix with dimension [num_points, iOut]
%          each row corresponds to the values of the hierarchical
%          surpluses associated with the basis function.
%
%       Note: The order and leading dimension must
%             match the points obtained form tsgGetPoints(...) command
%
% OUTPUT:
%
% The grid file associated with lGrid is modified
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC] = tsgMakeFilenames(lGrid.sName);

sCommand = [sTasGrid,' -loadhierarchy'];

sCommand = [sCommand, ' -gridfile ', sFileG];

tsgWriteMatrix(sFileV, mValues);
lClean.sFileV = 1;

sCommand = [sCommand, ' -vf ', sFileV];

[status, cmdout] = system(sCommand);

if (size(findstr('ERROR', cmdout)) ~= [0, 0])
    disp(cmdout);
    error('The tasgrid execurable returned an error, see above');
    return;
else
    if (~isempty(cmdout))
        fprintf(1,['Warning: Command had non-empty output:\n']);
        disp(cmdout);
    end
end

tsgCleanTempFiles(lGrid, lClean);

end
