function tsgSummary(lGrid)
%
% tsgSummary(lGrid)
%
% prints human readable summary of the grid properties
%
% INPUT:
%
% lGrid: a grid list created by tsgMakeXXX(...)
%
% OUTPUT:
%
% Text is printed out
% Note: lGrid is NOT modified by this command
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC] = tsgMakeFilenames(lGrid);

sCommand = [sTasGrid,' -summary'];

sCommand = [sCommand, ' -gridfile ', sFileG];

[status, cmdout] = system(sCommand);

if (max(size(strfind(cmdout, 'ERROR'))) ~= 0)
    disp(cmdout);
    error('The tasgrid execurable returned an error, see above');
    return;
else
    disp(cmdout);
end

end
