function [sFileG, sFileX, sFileV, sFileO, sFileW, sFileC, sFileL] = tsgMakeFilenames(lGrid)
%
% [sFileG, sFileX, sFileV, sFileO, sFileW, sFileC, sFileL] = tsgMakeFilenames(lGrid)
%
% given the paths from tsgGetPaths() it returns the filenames of all files associated with the name
%

[sFiles, sTasGrid] = tsgGetPaths();

if exist('getpid','builtin')
    id = getpid();
else
    id = feature('getpid');
end

if (isfield(lGrid, 'sFilename'))
    % new convention, the filename is stored into the lGrid object
    sFileG = lGrid.sFilename;
else
    % using the old naming convention
    sFileG = [sFiles, lGrid.sName,'_FileG']; % the filename to store the grid
end

sFileX = [sFiles, lGrid.sName,'.FileX',num2str(id)]; % file with points to evaluate the surrogate

sFileV = [sFiles, lGrid.sName,'.FileV']; % file with values to be loaded

sFileO = [sFiles, lGrid.sName,'.FileO',num2str(id)]; % file for output

sFileW = [sFiles, lGrid.sName,'.FileW']; % file with anisotropic weights

sFileC = [sFiles, lGrid.sName,'.FileC']; % file with the custom rule

sFileL = [sFiles, lGrid.sName,'.FileL']; % file with level limits

end
