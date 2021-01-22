function [sFileG, sFileX, sFileV, sFileO, sFileW, sFileC, sFileL] = tsgMakeFilenames(sGridName)
%
% [sFileG, sFileX, sFileV, sFileO, sFileW, sFileC, sFileL] = tsgMakeFilenames(sGridName)
%
% given the paths from tsgGetPaths() it returns the filenames of all files associated with the name
%

[sFiles, sTasGrid] = tsgGetPaths();

if exist('getpid','builtin')
    id = getpid();
else
    id = feature('getpid');
end

sFileG = [sFiles, sGridName,'_FileG']; % the filename to store the grid

sFileX = [sFiles, sGridName,'_FileX',num2str(id)]; % file with points to evaluate the surrogate

sFileV = [sFiles, sGridName,'_FileV']; % file with values to be loaded

sFileO = [sFiles, sGridName,'_FileO',num2str(id)]; % file for output

sFileW = [sFiles, sGridName,'_FileW']; % file with anisotropic weights

sFileC = [sFiles, sGridName,'_FileC']; % file with the custom rule

sFileL = [sFiles, sGridName,'_FileL']; % file with level limits

end
