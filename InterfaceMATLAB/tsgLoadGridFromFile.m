function [lGrid] = tsgLoadGridFromFile(sGridName, sFilename)
%
% [lGrid] = tsgLoadGridFromFile(sGridName, sFilename)
%
% Creates an lGrid object referencing an arbitrary file.
% The temporary folder will still be used for data such as points,
% weights and model values, therefore a unique name is still needed.
%
% The method is intended to simplify interfaces with other languages,
% e.g., the grid can be constructed with C++ and stored in a file
% then post-processing can done with MATLAB.
%
% Note: this does not check if there are already other lGrid objects
% associated with the file, aliasing can cause issues here.
%
% INPUT:
%
% sGridName: unique identifier string used for temporary filenames
%
% sFilename: existing file holding a Tasmanian sparse grid
%
% OUTPUT:
%
% lGrid: list containing information about the sparse grid,
%        can be used to call other functions
%

if (exist(sFilename, 'file') ~= 2)
    error(['invalid filename, cannot find file: ',sFilename]);
end

[sFiles, sTasGrid] = tsgGetPaths();
lGrid.sName = sGridName;
lGrid.sFilename = sFilename;

sCommand = [sTasGrid,' -summary -gridfile ', sFilename];

[status, cmdout] = system(sCommand);

if (max(size(strfind(cmdout, 'ERROR'))) ~= 0)
    disp(cmdout);
    error('The tasgrid execurable returned an error, see above');
end

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
