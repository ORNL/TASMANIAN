function lCustomRule = tsgReadCustomRuleFile(sFilename)
%
% reads a custom-tabulated rule from an ASCII file
%
% description: Test CT
% levels: 2
% 1 1
% 2 3
% 2.0 0.0
% 1.5 -0.5
% 0.5 0.5
%
% results in the list with properties
%
% sDescription = 'Test CT'
% iMaxLevel = 2
% vLevels = [1; 2]
% vPrecision = [1; 3]
% vNodes = [0.0, -0.5, 0.5]
% vWeights = [2.0, 1.5, 0.5]
%

sFilename = regexprep(sFilename, '\\ ', ' ');
fid = fopen(sFilename);

sLine = fgetl(fid);
if (strfind('description: ', sLine) ~= 1)
    error('Wrong file format of custom tabulated file on line 1!');
end
lCustomRule.sDescription = strrep(sLine, 'description: ', '');

sLine = fgetl(fid);
if (strfind('levels: ', sLine) ~= 1)
    error('Wrong file format of custom tabulated file on line 2!');
end
lCustomRule.iMaxLevel = str2double(strrep(sLine, 'levels: ', ''));

lCustomRule.vLevels = zeros(lCustomRule.iMaxLevel, 1);
lCustomRule.vPrecision = zeros(lCustomRule.iMaxLevel, 1);
for i=1:lCustomRule.iMaxLevel
    s = fscanf(fid, ' %f ', [1, 2]);
    lCustomRule.vLevels(i) = s(1);
    lCustomRule.vPrecision(i) = s(2);
end

lCustomRule.vNodes = zeros(sum(lCustomRule.vLevels), 1);
lCustomRule.vWeights = zeros(sum(lCustomRule.vLevels), 1);
for i=1:sum(lCustomRule.vLevels)
    s = fscanf(fid, ' %f ', [1, 2]);
    lCustomRule.vNodes(i) = s(2);
    lCustomRule.vWeights(i) = s(1);
end

fclose(fid);

end
