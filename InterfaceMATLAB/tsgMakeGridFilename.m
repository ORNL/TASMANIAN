function sGridFilename = tsgMakeGridFilename(sGridName)
%
% sGridFilename = tsgMakeGridFilename(sGridName)
%
% given the paths from tsgGetPaths() it returns the filename for the grid
%

[sFiles, sTasGrid] = tsgGetPaths();

if exist('getpid','builtin')
    id = getpid();
else
    id = feature('getpid');
end

sGridFilename = [sFiles, sGridName,'.tsgrid'];

end
