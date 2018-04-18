function [sFiles, sTasGrid] = tsgGetPaths()
%
% [sFiles, sTasGrid] = tsgGetPaths()
%
% You should edit the following variables accordingly
% sTasGrdid should be the path + executable of that tasgrid wrapper
% sFile should be root of the work files directory
%
% it is recommended to use absolute path

sTasGrid = ['@Tasmanian_tasgrid_path@'];
sFiles = ['@Tasmanian_MATLAB_WORK_FOLDER@'];

end
