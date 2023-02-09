function [ result ] = tsgUsingConstruction(lGrid)
%
% tsgUsingConstruction(lGrid)
%
% return 0/1 indicating whether dynamic construction is enabled
% corresponds to TasmanianSparseGrid::isUsingConstruction()
%
% INPUT:
%
% lGrid: a grid list created by tsgMakeXXX(...)
%
% OUTPUT:
%
% result will be either 1 if the dynamic construction has been enabled
%        or 0 if not enabled
% Note: lGrid is NOT modified by this command
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC] = tsgMakeFilenames(lGrid);

sCommand = [sTasGrid,' -using-construct'];

sCommand = [sCommand, ' -gridfile ', sFileG];

[status, cmdout] = system(sCommand);

if (max(size(strfind(cmdout, 'ERROR'))) ~= 0)
    disp(cmdout);
    error('The tasgrid execurable returned an error, see above');
    return;
else
    result = (norm( size(strfind(cmdout, 'enabled')) ) > 0);
end

end
