function tsgPlotPoints2D(points, fig)
%
% plotGrid2D(points, fig)
%
% plots two dimensional points in a figure
%
% INPUT:
%
% points: (matrix number_of_points X 2)
%          points created using tsgMakeXXX(...) command with iDim = 2
%
% fig: (optional integer figure number)
%      the plot will show on the given figure
%      by default everything is plotted on figure 1
%

if (nargin < 2)
    fig = 1;
end

if (size(points, 2) ~= 2)
    disp(['ERROR: size of q should be n by 2']);
	return
end

qm = min(points);
qM = max(points);

qD = (qM - qm);

qI = 0.05 * qD;

qm = qm - qI;
qM = qM + qI;

figure(fig)
hold on
for i = 1:size(points,1)
	pl = plot(points(i,1), points(i,2), '.');
    set(pl, 'MarkerSize', 15);
end

axis([qm(1),qM(1),qm(2),qM(2)]);

end