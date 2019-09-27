function tsgWriteCustomRuleFile(sFilename, lCustomRule)
%
% tsgWriteCustomRuleFile(sFilename, lCustomRule)
%
% uses the custom rule specified by the lCustomRule variable and writes
% it into a file with the given filename
%
% INPUT:
%
% sFilename: the file where the custom rule should be written
%
% lCustomRule:
%         the customRule structure that specifies the rule
%
% lCustomRule.sDescription:
%         the string describing the rule (for user information purposes only)
%
% lCustomRule.iMaxLevel:
%         the number (integer) of levels defined by the rule
%
% lCustomRule.vLevels:
%         a vector of length lCustomRule.iMaxLevel (integers)
%         defines the number of points on each level
%
% lCustomRule.vPrecision:
%         a vector of length lCustomRule.iMaxLevel (integers)
%         defines the precision (in polynomial order) or each level of the
%         quadrature rule, e.g., for Gauss-Legendre rule
%         custom_rule.precision = 2*custom_rule.levels-1
%         NOTE: grids of type other than ip* and qp* will ignore the
%               precision, in which case, precision can be set to a
%               vector of zeros
%
% lCustomRule.vNodes:
%         a vector of length sum(lCustomRule.vLevels) (floats)
%         for each level in order, defines the abscissas of the custom rule
%
% lCustomRule.vWeights:
%         a vector of length sum(lCustomRule.vLevels) (floats)
%         for each level in order, defines the weights of the custom rule
%         corresponding the lCustomRule.vNodes
%
% OUTPUT:
%         a file is created with the custom rule
%

sFilename = regexprep(sFilename, '\\ ', ' ');
fid = fopen(sFilename, 'w');

if (isfield(lCustomRule, 'sDescription'))
    sDesc = ['description: ', lCustomRule.sDescription,'\n'];
else
    sDesc = ['description: ', lCustomRule.description,'\n'];
end
fprintf(fid, sDesc);

if (isfield(lCustomRule, 'iMaxLevel'))
    sLevels = ['levels: ',num2str(lCustomRule.iMaxLevel),'\n'];
else
    sLevels = ['levels: ',num2str(lCustomRule.max_level),'\n'];
end

fprintf(fid, sLevels);

if (isfield(lCustomRule, 'vLevels'))
    for i = 1:lCustomRule.iMaxLevel
        fprintf(fid, '%d %d\n', [lCustomRule.vLevels(i), lCustomRule.vPrecision(i)]);
    end
else
    for i = 1:lCustomRule.max_level
        fprintf(fid, '%d %d\n', [lCustomRule.levels(i), lCustomRule.precision(i)]);
    end
end

if (isfield(lCustomRule, 'vWeights'))
    for i = 1:sum(lCustomRule.vLevels)
        fprintf(fid, '%2.20e %2.20e\n', [lCustomRule.vWeights(i), lCustomRule.vNodes(i)]);
    end
else
    for i = 1:sum(lCustomRule.levels)
        fprintf(fid, '%2.20e %2.20e\n', [lCustomRule.weights(i), lCustomRule.nodes(i)]);
    end
end

fclose(fid);

end
