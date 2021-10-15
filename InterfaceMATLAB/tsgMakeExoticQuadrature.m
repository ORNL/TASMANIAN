function lCustomRule = tsgMakeExoticQuadrature(iDepth, fShift, lWeightGrid, sDescription, bSymmetric)
%
% creates a custom-tabulated contain global exotic quadrature nodes and weights
%
% INPUT:
%
% iDepth: (integer non-negative)
%         controls the density of the grid
%
% fShift: (float)
%         the shift of the weight function
%
% lWeightGrid: (sparse grid of weight function)
%              list containing information about the surrogate/interpolant of the weight function; can be ASCII or binary format
%
% sDescription: (string)
%               description of the generated custom rule, e.g. "Exotic Quadrature 1" or "Sinc rule"
%
% bSymmetric: (bool)
%             if true, this declares that the weight function is symmetric
%
% OUTPUT:
%
% lCustomRule: (global grids of custom-tabulated rule)
%
%                structure defining the fields
%                   lCustomRule.sDescription
%                   lCustomRule.iMaxLevel
%                   lCustomRule.vLevels
%                   lCustomRule.vPrecision
%                   lCustomRule.vNodes
%                   lCustomRule.vWeights
%
%                  see help tsgWriteCustomRuleFile.m for definition of each field of the structure

if (isempty(sDescription))
    error('sDescription cannot be empty');
end

% generate filenames
[~, sTasGrid] = tsgGetPaths();
[~, ~, ~, ~, ~, sFileC, ~] = tsgMakeFilenames('');
[sWeightFileG, ~, ~, ~, ~, ~, ~] = tsgMakeFilenames(lWeightGrid.sName);

sCommand = [sTasGrid,' -makeexoquad'];
sCommand = [sCommand, ' -depth ',       num2str(iDepth)];
sCommand = [sCommand, ' -shift ',       num2str(fShift)];
sCommand = [sCommand, ' -weightfile ',  sWeightFileG];
sCommand = [sCommand, ' -description ', '"', sDescription, '"'];
if (bSymmetric)
    sCommand = [sCommand, ' -symmetric '];
end
sCommand = [sCommand, ' -outputfile ', sFileC];

[~, cmdout] = system(sCommand);

if (max(size(strfind('ERROR', cmdout))) ~= 0)
    disp(cmdout);
    error('The tasgrid execurable returned an error, see above');
end

% create lCustomRule object
lCustomRule = tsgReadCustomRuleFile(sFileC);

lClean.sFileC = 1;
lDummyGrid.sName = '';
tsgCleanTempFiles(lDummyGrid, lClean);

end
