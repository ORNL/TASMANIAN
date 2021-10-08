function lCustomRule = tsgReadCustomTabulated(sFilename)
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
lCustomRule.iMaxLevel = str2num(strrep(sLine, 'levels: ', ''));

lCustomRule.vLevels = NaN(lCustomRule.iMaxLevel, 1);
lCustomRule.vPrecision = NaN(lCustomRule.iMaxLevel, 1);
for i=1:lCustomRule.iMaxLevel
  sLine = fgetl(fid);
  aData = strsplit(sLine, ' ');
  lCustomRule.vLevels(i) = str2num(aData{1});
  lCustomRule.vPrecision(i) = str2num(aData{2});
end

lCustomRule.vNodes = NaN(sum(lCustomRule.vLevels), 1);
lCustomRule.vWeights = NaN(sum(lCustomRule.vLevels), 1);
for i=1:sum(lCustomRule.vLevels)
  sLine = fgetl(fid);
  aData = strsplit(sLine, ' ');
  lCustomRule.vNodes(i) = str2double(aData{2});
  lCustomRule.vWeights(i) = str2double(aData{1});
end

fclose(fid);

end
