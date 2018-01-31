function tsgMergeRefine(lGrid)
%
% tsgMergeRefine(lGrid)
%
% merges the points and loaded points and discards the loaded values
%
% INPUT:
%
% lGrid: a grid list created by tsgMake***(...)
%
% OUTPUT:
%
% The grid file associated with lGrid is modified
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC] = tsgMakeFilenames(lGrid.sName);

sCommand = [sTasGrid,' -mergerefine'];

sCommand = [sCommand, ' -gridfile ', sFileG];

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
end

end
