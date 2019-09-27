function [lGrid] = tsgReloadGrid(sGridName)
%
% [lGrid] = tsgReloadGrid(sGridName)
%
% Searches through the work files for files corresponding to a grid with
% name sName, then creates a new lGrid object associated with the grid.
%
% Note: this does not check if there are already other lGrid objects
% associated with the name, aliasing can cause issues here.
%
% INPUT:
%
% sGridName: name of a grid that exists in the work files
%
% OUTPUT:
%
% lGrid: list containing information about the sparse grid, can be used to
%        call other functions
%

[sFiles, sTasGrid] = tsgGetPaths();
[sFileG, sFileX, sFileV, sFileO, sFileW, sFileC] = tsgMakeFilenames(sGridName);

sTestFile = regexprep(sFileG, '\\ ', ' ');
if (exist(sTestFile, 'file') ~= 2)
    error(['There is no grid named: ',sGridName]);
end

sCommand = [sTasGrid,' -summary -gridfile ', sFileG];

[status, cmdout] = system(sCommand);

sLines = strsplit(cmdout, {'\n','\r'});

l1 = strsplit(sLines{2});

l2 = strsplit(sLines{3});

l3 = strsplit(sLines{4});

% create lGrid object
lGrid.sName = sGridName;

if (strcmp(l1{4}, 'Global'))
    lGrid.sType = 'global';
elseif (strcmp(l1{4}, 'Sequence'))
    lGrid.sType = 'sequence';
elseif (strcmp(l1{4}, 'Local'))
    lGrid.sType = 'localpolynomial';
elseif (strcmp(l1{4}, 'Wavelets'))
    lGrid.sType = 'wavelet';
elseif (strcmp(l1{4}, 'Fourier'))
    lGrid.sType = 'fourier';
end

lGrid.iDim = str2num(l2{3});
lGrid.iOut = str2num(l3{3});

end
