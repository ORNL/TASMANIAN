function tsgCoreTests()

% Indexes and polynomial spaces are experimental
% ListGridsByName(), Summary(), and PlotPoins2D() require human

disp(['']);
disp(['Testing TASMANIAN MATLAB interface']);
[sFiles, sTasGrid] = tsgGetPaths();
disp(['Tasmanian executable: ']);
disp(['  ',sTasGrid]);
disp(['Tasmanian work folder:']);
disp(['  ', sFiles]);
disp(['']);
[status, cmdout] = system([sTasGrid, ' -v']);
if (status ~= 0)
    disp(cmdout);
    error('There was an error while executing tasgrid.');
end
%k = 1;
%ll = 0;
%while(((k + 6) < length(cmdout)) && (ll < 9))
%    if ((cmdout(k) == ' ') && (cmdout(k+1) == ' ') && (cmdout(k+2) == ' ') && (cmdout(k+3) == ' '))
%        ll = ll + 1;
%        k = k + 4;
%    else
%        k = k + 1;
%    end
%end
%disp(cmdout(1:k));
disp(cmdout);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgMakeQuadrature()                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic quadrature calls and check for correct points and weight
[weights, points] = tsgMakeQuadrature(2, 'clenshaw-curtis', 'level', 1, 0);
tw = [4.0/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0]';
tp = [0.0, 0.0; 0.0, -1.0; 0.0, 1.0; -1.0, 0.0; 1.0, 0.0; ];
if ((norm(tw - weights) > 1.E-11) || (norm(tp - points) > 1.E-11))
    error('Mismatch in points and weights of simple quadrature 1');
end

[weights, points] = tsgMakeQuadrature(2, 'clenshaw-curtis', 'level', 2, 0);
if ((norm(points(4,2) + 1.0 / sqrt(2.0)) > 1.E-11) || (norm(points(5,2) - 1.0 / sqrt(2.0)) > 1.E-11) ...
    || (norm(weights(7) - 1.0 / 9.0) > 1.E-11))
    error('Mismatch in points and weights of simple quadrature 2');
end

[weights, points] = tsgMakeQuadrature(3, 'fejer2', 'level', 4, 0);
if ((norm(sum(weights) - 2.0^3) > 1.E-11) || (abs(sum(sum(points))) > 1.E-11))
    error('Mismatch in points and weights of simple quadrature 3');
end

[weights, points] = tsgMakeQuadrature(1, 'leja', 'level', 3, 0);
tw = [4.0/3.0, 1.0/3.0, 1.0/3.0, 0.0]';
tp = [0.0, 1.0, -1.0, sqrt(1.0/3.0)]';
if ((norm(tw - weights) > 1.E-11) || (norm(tp - points) > 1.E-11))
    error('Mismatch in points and weights of simple quadrature 4');
end

% test transform
[w, p] = tsgMakeQuadrature(3, 'clenshaw-curtis', 'level', 2, 0, [3.0 5.0; -7.0 -6.0; -12.0 17.0]);
if ((abs(norm(sum(w)) - 58.0) > 1.E-11) || (abs(max(p(:, 1)) - 5.0) > 1.E-11) ...
    || (abs(min(p(:, 1)) - 3.0) > 1.E-11) || (abs(max(p(:, 2)) + 6.0) > 1.E-11) ...
    || (abs(min(p(:, 2)) + 7.0) > 1.E-11) || (abs(max(p(:, 3)) - 17.0) > 1.E-11) ...
    || (abs(min(p(:, 3)) + 12.0) > 1.E-11))
    error('Mismatch in points and weights of simple quadrature: transform');
end

% test alpha/beta
[w, p] = tsgMakeQuadrature(1, 'gauss-hermite', 'level', 4, 0, [], [2.0;]);
if (abs(norm(sum(w)) - 0.5 * pi^0.5) > 1.E-11)
    error('Mismatch in points and weights of simple quadrature: alpha/beta');
end

% test anisotropy
[w, p] = tsgMakeQuadrature(2, 'leja', 'level', 2, 0, [], [], [2, 1]');
tp = [0.0 0.0; 0.0 1.0; 0.0 -1.0; 1.0 0.0;];
if ((abs(sum(w) - 4.0) > 1.E-11) || (norm(p - tp) > 1.E-11))
    error('Mismatch in points and weights of simple quadrature: anisotropy');
end

% test custom rule
vNodes = [];
vWeights = [];
for iL = 0:4
    [w, p] = tsgMakeQuadrature(1, 'gauss-legendre', 'level', iL, -1);
    vNodes = [vNodes; p];
    vWeights = [vWeights; w];
end
lCustomRule.sDescription = 'Test Gauss-Legendre';
lCustomRule.iMaxLevel = 5;
lCustomRule.vLevels = 1:5;
lCustomRule.vPrecision = 2 * lCustomRule.vLevels - 1;
lCustomRule.vNodes = vNodes;
lCustomRule.vWeights = vWeights;
[tw, tp] = tsgMakeQuadrature(2, 'gauss-legendre', 'level', 3, 0, [], [], [2, 1]');
[w, p] = tsgMakeQuadrature(2, 'custom-tabulated', 'level', 3, 0, [], [], [2, 1]', lCustomRule);
if ((norm(tw - w) > 1.E-11) || (norm(tp - p) > 1.E-11))
    error('Mismatch in points and weights of simple quadrature: custom rule');
end

% test conformal
for iL = 7:8
    [w, p] = tsgMakeQuadrature(2, 'gauss-patterson', 'qptotal', iL, -1);
    [wc, pc] = tsgMakeQuadrature(2, 'gauss-patterson', 'qptotal', iL, -1, [], [], [], [], 'asin', [4, 4]);

    I  = sum(w  .* (1.0 ./ ((1.0 + 5.0 .* p(:,1).^2)  .* (1.0 + 5.0 .* p(:,2).^2))));
    Ic = sum(wc .* (1.0 ./ ((1.0 + 5.0 .* pc(:,1).^2) .* (1.0 + 5.0 .* pc(:,2).^2))));
    %[abs(I - 1.028825601981092^2), abs(Ic - 1.028825601981092^2)]

    if (abs(I - 1.028825601981092^2) < abs(Ic - 1.028825601981092^2))
        error('Mismatch in points and weights of simple quadrature: conformal map');
    end
end

% test level limits
[w, p] = tsgMakeQuadrature(2, 'clenshaw-curtis', 'qptotal', 20, -1, [], [], [], [], [], [], [1, 3]);
[tw, tp] = tsgMakeQuadrature(2, 'clenshaw-curtis', 'tensor', 1, -1, [], [], [1, 3], [], [], [], []);
% when the limitation of iDepth far exceeds the level limits, the grid converges to a full tensor grid
if (abs(sum(w) - 4.0) > 1.E-11)
    error('Mismatch in points and weights of simple quadrature: level limit, sum of weights');
end
if ((norm(w - tw) > 1.E-11) || (norm(p - tp) > 1.E-11))
    error('Mismatch in points and weights of simple quadrature: level limit, points and weights');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 tsgMakeExoticQuadrature()                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a sinc function surrogate.
[lWeightGrid, points] = tsgMakeGlobal('Sinc_Surrogate', 1, 1, 'gauss-legendre', 'qptotal', 200);
vals = zeros(length(points), 1);
for i=1:length(points)
    if abs(points(i)) <= 1e-16
        vals(i) = 1.0;
    else
        vals(i) = sin(points(i)) / points(i);
    end
end
tsgLoadValues(lWeightGrid, vals);

% Create the exotic quadrature rule and test it on 1D and 2D instances.
iDepth = 40;
lCustomRule = tsgMakeExoticQuadrature(iDepth, 1.0, lWeightGrid, 'Sinc_Exoquad', true);
[w1, p1] = tsgMakeQuadrature(1, 'custom-tabulated', 'qptotal', iDepth, 0, [], [], [], lCustomRule);
I1 = sum(w1 .* exp(-p1 .* p1));
if (abs(I1 - 1.4321357541271255) > 1E-11)
    error('Mismatch in generated sinc-weighted integral of tsgMakeExoticQuadrature() for dimension 1');
end
[w2, p2] = tsgMakeQuadrature(2, 'custom-tabulated', 'qptotal', iDepth, 0, [], [], [], lCustomRule);
I2 = sum(w2 .* exp(-p2(:, 1) .* p2(:, 1) - p2(:, 2) .* p2(:, 2)));
if (abs(I2 - 1.4321357541271255 ^ 2) > 1E-11)
    error('Mismatch in generated sinc-weighted integral of tsgMakeExoticQuadrature() for dimension 2');
end

% Clean up
tsgDeleteGrid(lWeightGrid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgMakeGlobal()                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lGrid, p] = tsgMakeGlobal('_tsgcoretests_lgrid', 2, 1, 'clenshaw-curtis', 'level', 1);
tp = [0.0, 0.0; 0.0, -1.0; 0.0, 1.0; -1.0, 0.0; 1.0, 0.0;];
if (norm(tp - p) > 1.E-11)
    error('Mismatch in tsgMakeGlobal: core case 1');
end

% test transform
[lGrid, p] = tsgMakeGlobal('_tsgcoretests_lgrid', 3, 1, 'clenshaw-curtis', 'level', 2, [3.0 5.0; -7.0 -6.0; -12.0 17.0]);
if ((abs(max(p(:, 1)) - 5.0) > 1.E-11) || (abs(min(p(:, 1)) - 3.0) > 1.E-11) ...
    || (abs(max(p(:, 2)) + 6.0) > 1.E-11) || (abs(min(p(:, 2)) + 7.0) > 1.E-11) ...
    || (abs(max(p(:, 3)) - 17.0) > 1.E-11) || (abs(min(p(:, 3)) + 12.0) > 1.E-11))
    error('Mismatch in tsgMakeGlobal: transform');
end
[w, p] = tsgGetQuadrature(lGrid);
if ((abs(norm(sum(w)) - 58.0) > 1.E-11) || (abs(max(p(:, 1)) - 5.0) > 1.E-11) ...
    || (abs(min(p(:, 1)) - 3.0) > 1.E-11) || (abs(max(p(:, 2)) + 6.0) > 1.E-11) ...
    || (abs(min(p(:, 2)) + 7.0) > 1.E-11) || (abs(max(p(:, 3)) - 17.0) > 1.E-11) ...
    || (abs(min(p(:, 3)) + 12.0) > 1.E-11))
    error('Mismatch in tsgMakeGlobal: getQuadrature');
end

% test alpha/beta
[lGrid] = tsgMakeGlobal('_tsgcoretests_lgrid', 1, 1, 'gauss-hermite', 'level', 4, [], [2.0;]);
[w, p] = tsgGetQuadrature(lGrid);
if (abs(norm(sum(w)) - 0.5 * pi^0.5) > 1.E-11)
    error('Mismatch in tsgMakeGlobal: alpha/beta');
end

% test anisotropy
[lGrid] = tsgMakeGlobal('_tsgcoretests_lgrid', 2, 1, 'leja', 'level', 2, [], [], [2, 1]);
w = []; p = [];
[w, p] = tsgGetQuadrature(lGrid);
tp = [0.0 0.0; 0.0 1.0; 0.0 -1.0; 1.0 0.0;];
if ((abs(sum(w) - 4.0) > 1.E-11) || (norm(p - tp) > 1.E-11))
    error('Mismatch in tsgMakeGlobal: anisotropy');
end
[lGrid] = tsgMakeGlobal('_tsgcoretests_lgrid', 2, 1, 'leja', 'level', 2, [], [], [2, 1]');
w = []; p = [];
[w, p] = tsgGetQuadrature(lGrid);
tp = [0.0 0.0; 0.0 1.0; 0.0 -1.0; 1.0 0.0;];
if ((abs(sum(w) - 4.0) > 1.E-11) || (norm(p - tp) > 1.E-11))
    error('Mismatch in tsgMakeGlobal: anisotropy');
end

% test custom rule
vNodes = [];
vWeights = [];
for iL = 0:4
    [w, p] = tsgMakeQuadrature(1, 'gauss-legendre', 'level', iL, -1);
    vNodes = [vNodes; p];
    vWeights = [vWeights; w];
end
lCustomRule.sDescription = 'Test Gauss-Legendre';
lCustomRule.iMaxLevel = 5;
lCustomRule.vLevels = 1:5;
lCustomRule.vPrecision = 2 * lCustomRule.vLevels - 1;
lCustomRule.vNodes = vNodes;
lCustomRule.vWeights = vWeights;
[tw, tp] = tsgMakeQuadrature(2, 'gauss-legendre', 'level', 3, 0, [], [], [2, 1]');
[lGrid] = tsgMakeGlobal('_tsgcoretests_lgrid', 2, 1, 'custom-tabulated', 'level', 3, [], [], [2, 1]', lCustomRule);
w = []; p = [];
[w, p] = tsgGetQuadrature(lGrid);
if ((norm(tw - w) > 1.E-11) || (norm(tp - p) > 1.E-11))
    error('Mismatch in tsgMakeGlobal: custom rule');
end

% test conformal
pnts = [-1.0 + 2.0 * rand(100, 2)];
tres = 1.0 ./ ((1.0 + 5.0 .* pnts(:,1).^2) .* (1.0 + 5.0 .* pnts(:,2).^2));
for iL = 7:8
    [lGrid, p] = tsgMakeGlobal('_tsgcoretests_lgrid', 2, 1, 'clenshaw-curtis', 'iptotal', iL, [], [], [], [], 'asin', [4, 4]);
    n1 = size(p, 1);
    v = 1.0 ./ ((1.0 + 5.0 .* p(:,1).^2) .* (1.0 + 5.0 .* p(:,2).^2));
    tsgLoadValues(lGrid, v);
    [resc] = tsgEvaluate(lGrid, pnts);
    p = [];
    [lGrid, p] = tsgMakeGlobal('_tsgcoretests_lgrid', 2, 1, 'clenshaw-curtis', 'iptotal', iL);
    v = 1.0 ./ ((1.0 + 5.0 .* p(:,1).^2) .* (1.0 + 5.0 .* p(:,2).^2));
    tsgLoadValues(lGrid, v);
    [resn] = tsgEvaluate(lGrid, pnts);
    if (size(p, 1) ~= n1)
        error('Mismatch in tsgMakeGlobal: conformal number of points');
    end
    if (norm(tres - resc) > norm(tres - resn))
        error('Mismatch in tsgMakeGlobal: conformal error');
    end
end

% test level limits
[lGrid, p1] = tsgMakeGlobal('_tsgcoretests_lgrid', 2, 1, 'clenshaw-curtis', 'iptotal', 20, [], [], [], [], [], [], [1, 3]);
[lGrid, p2] = tsgMakeGlobal('_tsgcoretests_lgrid', 2, 1, 'clenshaw-curtis', 'tensor', 1, [], [], [1, 3]);
if (norm(p1 - p2) > 1.E-11)
    error('Mismatch in tsgMakeGlobal: level limits');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgDeleteGrid()                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sFiles, sTasGrid] = tsgGetPaths();
sFiles = regexprep(sFiles, '\\ ', ' ');
if (~exist([sFiles,'_tsgcoretests_lgrid_FileG'], 'file'))
    error('Mismatch in tsgDeleteGrid: cannot find a file that should exist');
end
tsgDeleteGrid(lGrid);
if (exist([sFiles,'_tsgcoretests_lgrid_FileG'], 'file'))
    error('Mismatch in tsgDeleteGrid: did not delete the file');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgCopyGrid()                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lGridA, p] = tsgMakeGlobal('_tsgcoretests_lgridA', 2, 1, 'clenshaw-curtis', 'iptotal', 3, [], [], [], [], 'asin', [4, 4]);
[lGridB] = tsgCopyGrid(lGridA, '_tsgcoretests_lgridB');
tsgDeleteGrid(lGridA);
[p2] = tsgGetPoints(lGridB);
if (norm(p - p2) > 1.E-11)
    error('Mismatch in tsgCopyGrid: did not delete the file');
end
tsgDeleteGrid(lGridB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgMakeSequence()                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lGrid, p] = tsgMakeSequence('_tsgcoretests_lgrid2', 2, 1, 'min-lebesgue', 'level', 3);
tp = [0.0, 0.0; 0.0, 1.0; 0.0, -1.0; 0.0, sqrt(1.0/3.0); 1.0, 0.0; 1.0, 1.0; 1.0, -1.0; -1.0, 0.0; -1.0, 1.0; sqrt(1.0/3.0), 0.0;];
if (norm(tp - p) > 1.E-11)
    error('Mismatch in tsgMakeSequence: core case 1');
end

% test transform
[lGrid, p] = tsgMakeSequence('_tsgcoretests_lgrid2', 3, 1, 'rleja', 'level', 2, [3.0 5.0; -7.0 -6.0; -12.0 17.0]);
if ((abs(max(p(:, 1)) - 5.0) > 1.E-11) || (abs(min(p(:, 1)) - 3.0) > 1.E-11) ...
    || (abs(max(p(:, 2)) + 6.0) > 1.E-11) || (abs(min(p(:, 2)) + 7.0) > 1.E-11) ...
    || (abs(max(p(:, 3)) - 17.0) > 1.E-11) || (abs(min(p(:, 3)) + 12.0) > 1.E-11))
    error('Mismatch in tsgMakeSequence: transform');
end
[w, p] = tsgGetQuadrature(lGrid);
if ((abs(norm(sum(w)) - 58.0) > 1.E-11) || (abs(max(p(:, 1)) - 5.0) > 1.E-11) ...
    || (abs(min(p(:, 1)) - 3.0) > 1.E-11) || (abs(max(p(:, 2)) + 6.0) > 1.E-11) ...
    || (abs(min(p(:, 2)) + 7.0) > 1.E-11) || (abs(max(p(:, 3)) - 17.0) > 1.E-11) ...
    || (abs(min(p(:, 3)) + 12.0) > 1.E-11))
    error('Mismatch in tsgMakeSequence: getQuadrature');
end

% test anisotropy
[lGrid] = tsgMakeSequence('_tsgcoretests_lgrid2', 2, 1, 'leja', 'level', 2, [], [2, 1]);
w = []; p = [];
[w, p] = tsgGetQuadrature(lGrid);
tp = [0.0 0.0; 0.0 1.0; 0.0 -1.0; 1.0 0.0;];
if ((abs(sum(w) - 4.0) > 1.E-11) || (norm(p - tp) > 1.E-11))
    error('Mismatch in tsgMakeSequence: anisotropy');
end
[lGrid] = tsgMakeSequence('_tsgcoretests_lgrid2', 2, 1, 'leja', 'level', 2, [], [2, 1]');
w = []; p = [];
[w, p] = tsgGetQuadrature(lGrid);
tp = [0.0 0.0; 0.0 1.0; 0.0 -1.0; 1.0 0.0;];
if ((abs(sum(w) - 4.0) > 1.E-11) || (norm(p - tp) > 1.E-11))
    error('Mismatch in tsgMakeSequence: anisotropy');
end

% test conformal
pnts = [-1.0 + 2.0 * rand(100, 2)];
tres = 1.0 ./ ((1.0 + 5.0 .* pnts(:,1).^2) .* (1.0 + 5.0 .* pnts(:,2).^2));
for iL = 7:8
    [lGrid, p] = tsgMakeSequence('_tsgcoretests_lgrid', 2, 1, 'rleja', 'iptotal', iL, [], [], 'asin', [4, 4]);
    n1 = size(p, 1);
    v = 1.0 ./ ((1.0 + 5.0 .* p(:,1).^2) .* (1.0 + 5.0 .* p(:,2).^2));
    tsgLoadValues(lGrid, v);
    [resc] = tsgEvaluate(lGrid, pnts);
    p = [];
    [lGrid, p] = tsgMakeSequence('_tsgcoretests_lgrid', 2, 1, 'rleja', 'iptotal', iL);
    v = 1.0 ./ ((1.0 + 5.0 .* p(:,1).^2) .* (1.0 + 5.0 .* p(:,2).^2));
    tsgLoadValues(lGrid, v);
    [resn] = tsgEvaluate(lGrid, pnts);
    if (size(p, 1) ~= n1)
        error('Mismatch in tsgMakeSequence: conformal number of points');
    end
    if (norm(tres - resc) > norm(tres - resn))
        error('Mismatch in tsgMakeSequence: conformal error');
    end
end

% test level limits
[lGrid, p1] = tsgMakeSequence('_tsgcoretests_lgrid', 2, 1, 'min-delta', 'iptotal', 20, [], [], [], [], [1, 3]);
[lGrid, p2] = tsgMakeSequence('_tsgcoretests_lgrid', 2, 1, 'min-delta', 'tensor', 1, [], [1, 3]);
if (norm(p1 - p2) > 1.E-11)
    error('Mismatch in tsgMakeSequence: level limits');
end
tsgDeleteGrid(lGrid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgDeleteGridByName()                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sFiles, sTasGrid] = tsgGetPaths();
sFiles = regexprep(sFiles, '\\ ', ' ');
if (~exist([sFiles,'_tsgcoretests_lgrid2_FileG'], 'file'))
    error('Mismatch in tsgDeleteGrid: cannot find a file that should exist');
end
tsgDeleteGridByName('_tsgcoretests_lgrid2');
if (exist([sFiles,'_tsgcoretests_lgrid2_FileG'], 'file'))
    error('Mismatch in tsgDeleteGrid: did not delete the file');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgMakeLocalPolynomial()                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test transform
[lGrid, p] = tsgMakeLocalPolynomial('_tsgcoretests_lgrid2', 3, 1, 'localp', 2, 2, [3.0 5.0; -7.0 -6.0; -12.0 17.0]);
if ((abs(max(p(:, 1)) - 5.0) > 1.E-11) || (abs(min(p(:, 1)) - 3.0) > 1.E-11) ...
    || (abs(max(p(:, 2)) + 6.0) > 1.E-11) || (abs(min(p(:, 2)) + 7.0) > 1.E-11) ...
    || (abs(max(p(:, 3)) - 17.0) > 1.E-11) || (abs(min(p(:, 3)) + 12.0) > 1.E-11))
    error('Mismatch in tsgMakeLocalPolynomial: transform');
end
[w, p] = tsgGetQuadrature(lGrid);
if ((abs(norm(sum(w)) - 58.0) > 1.E-11) || (abs(max(p(:, 1)) - 5.0) > 1.E-11) ...
    || (abs(min(p(:, 1)) - 3.0) > 1.E-11) || (abs(max(p(:, 2)) + 6.0) > 1.E-11) ...
    || (abs(min(p(:, 2)) + 7.0) > 1.E-11) || (abs(max(p(:, 3)) - 17.0) > 1.E-11) ...
    || (abs(min(p(:, 3)) + 12.0) > 1.E-11))
    error('Mismatch in tsgMakeLocalPolynomial: getQuadrature');
end
tsgDeleteGrid(lGrid);

% test conformal
pnts = [-1.0 + 2.0 * rand(100, 2)];
tres = 1.0 ./ ((1.0 + 5.0 .* pnts(:,1).^2) .* (1.0 + 5.0 .* pnts(:,2).^2));
for iL = 3:4
    [lGrid, p] = tsgMakeLocalPolynomial('_tsgcoretests_lgrid', 2, 1, 'semi-localp', iL, 2, [], 'asin', [4, 4]);
    n1 = size(p, 1);
    v = 1.0 ./ ((1.0 + 5.0 .* p(:,1).^2) .* (1.0 + 5.0 .* p(:,2).^2));
    tsgLoadValues(lGrid, v);
    [resc] = tsgEvaluate(lGrid, pnts);
    p = [];
    [lGrid, p] = tsgMakeLocalPolynomial('_tsgcoretests_lgrid', 2, 1, 'semi-localp', iL, 2);
    v = 1.0 ./ ((1.0 + 5.0 .* p(:,1).^2) .* (1.0 + 5.0 .* p(:,2).^2));
    tsgLoadValues(lGrid, v);
    [resn] = tsgEvaluate(lGrid, pnts);
    if (size(p, 1) ~= n1)
        error('Mismatch in tsgMakeLocalPolynomial: conformal number of points');
    end
    if (norm(tres - resc) > norm(tres - resn))
        error('Mismatch in tsgMakeLocalPolynomial: conformal error');
    end
end
tsgDeleteGrid(lGrid);
% polynomial order is tested in tsgEvaluate()

% level limits
[lGrid, p] = tsgMakeLocalPolynomial('_tsgcoretests_lgrid', 3, 1, 'semi-localp', 3, 2, [], [], [], [1, 2, 3]);
if (min(abs(p(:,1) - 0.5)) < 1.E-8)
    error('Mismatch in tsgMakeLocalPolynomial: level limit, dim 1');
end
if (min(abs(p(:,2) - 0.75)) < 1.E-8)
    error('Mismatch in tsgMakeLocalPolynomial: level limit, dim 2');
end
if (min(abs(p(:,3) - 0.125)) < 1.E-8)
    error('Mismatch in tsgMakeLocalPolynomial: level limit, dim 3');
end
tsgDeleteGrid(lGrid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgMakeWavelet()                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test transform
[lGrid, p] = tsgMakeWavelet('_tsgcoretests_lgrid', 3, 1, 2, 1, [3.0 5.0; -7.0 -6.0; -12.0 17.0]);
if ((abs(max(p(:, 1)) - 5.0) > 1.E-11) || (abs(min(p(:, 1)) - 3.0) > 1.E-11) ...
    || (abs(max(p(:, 2)) + 6.0) > 1.E-11) || (abs(min(p(:, 2)) + 7.0) > 1.E-11) ...
    || (abs(max(p(:, 3)) - 17.0) > 1.E-11) || (abs(min(p(:, 3)) + 12.0) > 1.E-11))
    error('Mismatch in tsgMakeWavelet: transform');
end
[w, p] = tsgGetQuadrature(lGrid);
if ((abs(norm(sum(w)) - 58.0) > 1.E-11) || (abs(max(p(:, 1)) - 5.0) > 1.E-11) ...
    || (abs(min(p(:, 1)) - 3.0) > 1.E-11) || (abs(max(p(:, 2)) + 6.0) > 1.E-11) ...
    || (abs(min(p(:, 2)) + 7.0) > 1.E-11) || (abs(max(p(:, 3)) - 17.0) > 1.E-11) ...
    || (abs(min(p(:, 3)) + 12.0) > 1.E-11))
    error('Mismatch in tsgMakeWavelet: getQuadrature');
end

% correctness of 1-D
[lGrid, pw] = tsgMakeWavelet('_tsgcoretests_lgrid', 1, 1, 2, 1);
[lGrid, p] = tsgMakeLocalPolynomial('_tsgcoretests_lgrid', 1, 1, 'localp', 3, 1);
if (norm(pw - p) > 1.E-11)
    error('Mismatch in tsgMakeWavelet: points');
end

% test conformal
for iL = 3:4
    [lGrid, p] = tsgMakeWavelet('_tsgcoretests_lgrid', 2, 1, iL, 1, [], 'asin', [4, 4]);
    [w, p] = tsgGetQuadrature(lGrid);
    [lGrid, p] = tsgMakeWavelet('_tsgcoretests_lgrid', 2, 1, iL, 1);
    [wc, pc] = tsgGetQuadrature(lGrid);

    I  = sum(w  .* (1.0 ./ ((1.0 + 5.0 .* p(:,1).^2)  .* (1.0 + 5.0 .* p(:,2).^2))));
    Ic = sum(wc .* (1.0 ./ ((1.0 + 5.0 .* pc(:,1).^2) .* (1.0 + 5.0 .* pc(:,2).^2))));

    if (abs(I - 1.028825601981092^2) < abs(Ic - 1.028825601981092^2))
        error('Mismatch in points and weights of simple quadrature: conformal map');
    end

end
tsgDeleteGrid(lGrid);

% level limits
[lGrid, p] = tsgMakeWavelet('_tsgcoretests_lgrid', 3, 1, 2, 1, [], [], [], [0, 1, 2]);
if (min(abs(p(:,1) - 0.5)) < 1.E-8)
    error('Mismatch in tsgMakeLocalPolynomial: level limit, dim 1');
end
if (min(abs(p(:,2) - 0.75)) < 1.E-8)
    error('Mismatch in tsgMakeLocalPolynomial: level limit, dim 2');
end
if (min(abs(p(:,3) - 0.125)) < 1.E-8)
    error('Mismatch in tsgMakeLocalPolynomial: level limit, dim 3');
end
tsgDeleteGrid(lGrid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgMakeFourier()                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test transform
[lGrid, p] = tsgMakeFourier('_tsgcoretests_lgrid', 3, 1, 'level', 1, [3.0 5.0; -7.0 -6.0; -12.0 17.0]);
if ((abs(max(p(:, 1)) - 13.0/3.0) > 1.E-11) || (abs(min(p(:, 1)) - 3.0) > 1.E-11) ...
    || (abs(max(p(:, 2)) + 19.0/3.0) > 1.E-11) || (abs(min(p(:, 2)) + 7.0) > 1.E-11) ...
    || (abs(max(p(:, 3)) - 22.0/3.0) > 1.E-11) || (abs(min(p(:, 3)) + 12.0) > 1.E-11))
    error('Mismatch in tsgMakeFourier: transform');
end
[w, p] = tsgGetQuadrature(lGrid);
if ((abs(norm(sum(w)) - 58.0) > 1.E-11) || (abs(max(p(:, 1)) - 13.0/3.0) > 1.E-11) ...
    || (abs(min(p(:, 1)) - 3.0) > 1.E-11) || (abs(max(p(:, 2)) + 19.0/3.0) > 1.E-11) ...
    || (abs(min(p(:, 2)) + 7.0) > 1.E-11) || (abs(max(p(:, 3)) - 22.0/3.0) > 1.E-11) ...
    || (abs(min(p(:, 3)) + 12.0) > 1.E-11))
    error('Mismatch in tsgMakeFourier: getQuadrature');
end

% correctness of 1-D
[lGrid, p] = tsgMakeFourier('_tsgcoretests_lgrid', 1, 1, 'level', 2);
tp = [0.0; 1.0/3.0; 2.0/3.0; 1.0/9.0; 2.0/9.0; 4.0/9.0; 5.0/9.0; 7.0/9.0; 8.0/9.0];
if (norm(p - tp) > 1.E-11)
    error('Mismatch in tsgMakeFourier: points');
end

% level limits
[lGrid, p] = tsgMakeFourier('_tsgcoretests_lgrid', 3, 1, 'level', 3, [], [], [], [], [0, 1, 2]);
if (max(abs(p(:,1))) > 1.E-8)
    error('Mismatch in tsgMakeFourier: level limit, dim 1');
end
if (min(abs(p(:,2) - 1.0/9.0)) < 1.E-8)
    error('Mismatch in tsgMakeFourier: level limit, dim 2');
end
if (min(abs(p(:,3) - 1.0/27.0)) < 1.E-8)
    error('Mismatch in tsgMakeFourier: level limit, dim 3');
end
tsgDeleteGrid(lGrid);

disp(['tsgMake* functions:       PASS']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgGetPoints()                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lGrid] = tsgMakeGlobal('_tsgcoretests_lgrid', 2, 1, 'clenshaw-curtis', 'level', 1);
tp = [0.0, 0.0; 0.0, -1.0; 0.0, 1.0; -1.0, 0.0; 1.0, 0.0;];
p = [];
[p] = tsgGetPoints(lGrid);
if (norm(tp - p) > 1.E-11)
    error('Mismatch in tsgGetPoints: core case 1');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgGetNeededPoints()                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% covered in tsgLoadValues() and tsgRefine*()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgGetQuadrature()                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% covered in tsgMakeGlobal() and tsgMakeSequence()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgLoadValues()                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lGrid, p] = tsgMakeGlobal('_tsgcoretests_lgrid', 2, 2, 'min-delta', 'level', 4);
[pn] = tsgGetNeededPoints(lGrid);
if (norm(p - pn) > 1.E-11)
    error('Mismatch in tsgLoadValues: tsgGetNeededPoints case 1');
end
v = [exp(-p(:,1).^2 -p(:,2).^2), cos(-p(:,1) -2.0 * p(:,2))];
tsgLoadValues(lGrid, v);
[pn] = tsgGetPoints(lGrid);
if (norm(p - pn) > 1.E-11)
    error('Mismatch in tsgLoadValues: tsgGetPoints');
end

[pn] = tsgGetNeededPoints(lGrid);
if (max(size(pn)) ~= 0)
    error('Mismatch in tsgLoadValues: tsgGetNeededPoints case 2');
end
tsgDeleteGrid(lGrid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgEvaluate()                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lGrid, p] = tsgMakeLocalPolynomial('_tsgcoretests_lp', 2, 4, 'localp', 1, 1);
v = [0.3*ones(size(p,1), 1), p(:,1) + p(:,2), ...
     p(:,1).^2 + p(:,2).^2 + p(:,1).*p(:,2), ...
     p(:,1).^3 + p(:,2).^3 + p(:,1).*(p(:,2).^2) ];
tsgLoadValues(lGrid, v);
p = [1.0/3.0, 1.0/3.0; pi/6.0, -sqrt(2.0)/2.0;];
tv = [0.3*ones(size(p,1), 1), p(:,1) + p(:,2), ...
      p(:,1).^2 + p(:,2).^2 + p(:,1).*p(:,2), ...
      p(:,1).^3 + p(:,2).^3 + p(:,1).*(p(:,2).^2) ];
[res] = tsgEvaluate(lGrid, p);
for i = 1:2
    if (norm(res(:,i) - tv(:,i)) > 1.E-11)
        error(['Mismatch in tsgEvaluate: case 1, output ',num2str(i)]);
    end
end
for i = 3:4
    if (norm(res(:,i) - tv(:,i)) < 1.E-8)
        error(['Mismatch in tsgEvaluate: case 1, output ',num2str(i)]);
    end
end

[lGrid, p] = tsgMakeLocalPolynomial('_tsgcoretests_lp', 2, 4, 'localp', 1, 2);
v = [0.3*ones(size(p,1), 1), p(:,1) + p(:,2), ...
     p(:,1).^2 + p(:,2).^2 + p(:,1).*p(:,2), ...
     p(:,1).^3 + p(:,2).^3 + p(:,1).*(p(:,2).^2) ];
tsgLoadValues(lGrid, v);
p = [1.0/3.0, 1.0/3.0; pi/6.0, -sqrt(2.0)/2.0;];
tv = [0.3*ones(size(p,1), 1), p(:,1) + p(:,2), ...
      p(:,1).^2 + p(:,2).^2 + p(:,1).*p(:,2), ...
      p(:,1).^3 + p(:,2).^3 + p(:,1).*(p(:,2).^2) ];
[res] = tsgEvaluate(lGrid, p);
for i = 1:2
    if (norm(res(:,i) - tv(:,i)) > 1.E-11)
        error(['Mismatch in tsgEvaluate: case 2, output ',num2str(i)]);
    end
end
for i = 3:4
    if (norm(res(:,i) - tv(:,i)) < 1.E-8)
        error(['Mismatch in tsgEvaluate: case 2, output ',num2str(i)]);
    end
end

[lGrid, p] = tsgMakeLocalPolynomial('_tsgcoretests_lp', 2, 4, 'semi-localp', 1, 2);
v = [0.3*ones(size(p,1), 1), p(:,1) + p(:,2), ...
     p(:,1).^2 + p(:,2).^2, ...
     p(:,1).^3 + p(:,2).^3 + p(:,1).*(p(:,2).^2) ];
tsgLoadValues(lGrid, v);
p = [1.0/3.0, 1.0/3.0; pi/6.0, -sqrt(2.0)/2.0;];
tv = [0.3*ones(size(p,1), 1), p(:,1) + p(:,2), ...
      p(:,1).^2 + p(:,2).^2, ...
      p(:,1).^3 + p(:,2).^3 + p(:,1).*(p(:,2).^2) ];
[res] = tsgEvaluate(lGrid, p);
for i = 1:3
    if (norm(res(:,i) - tv(:,i)) > 1.E-11)
        error(['Mismatch in tsgEvaluate: case 3, output ',num2str(i)]);
    end
end
for i = 4:4
    if (norm(res(:,i) - tv(:,i)) < 1.E-8)
        error(['Mismatch in tsgEvaluate: case 3, output ',num2str(i)]);
    end
end

[lGrid, p] = tsgMakeLocalPolynomial('_tsgcoretests_lp', 2, 4, 'localp', 2, 2);
v = [0.3*ones(size(p,1), 1), p(:,1) + p(:,2), ...
     p(:,1).^2 + p(:,2).^2 + p(:,1).*p(:,2), ...
     p(:,1).^3 + p(:,2).^3 + p(:,1).*(p(:,2).^2) ];
tsgLoadValues(lGrid, v);
p = [1.0/3.0, 1.0/3.0; pi/6.0, -sqrt(2.0)/2.0;];
tv = [0.3*ones(size(p,1), 1), p(:,1) + p(:,2), ...
      p(:,1).^2 + p(:,2).^2 + p(:,1).*p(:,2), ...
      p(:,1).^3 + p(:,2).^3 + p(:,1).*(p(:,2).^2) ];
[res] = tsgEvaluate(lGrid, p);
for i = 1:3
    if (norm(res(:,i) - tv(:,i)) > 1.E-11)
        error(['Mismatch in tsgEvaluate: case 4, output ',num2str(i)]);
    end
end
for i = 4:4
    if (norm(res(:,i) - tv(:,i)) < 1.E-8)
        error(['Mismatch in tsgEvaluate: case 4, output ',num2str(i)]);
    end
end

[lGrid, p] = tsgMakeLocalPolynomial('_tsgcoretests_lp', 2, 4, 'localp', 3, 3);
v = [0.3*ones(size(p,1), 1), p(:,1) + p(:,2), ...
     p(:,1).^2 + p(:,2).^2 + p(:,1).*p(:,2), ...
     p(:,1).^3 + p(:,2).^3 + p(:,1).*(p(:,2).^2) ];
tsgLoadValues(lGrid, v);
p = [1.0/3.0, 1.0/3.0; pi/6.0, -sqrt(2.0)/2.0;];
tv = [0.3*ones(size(p,1), 1), p(:,1) + p(:,2), ...
      p(:,1).^2 + p(:,2).^2 + p(:,1).*p(:,2), ...
      p(:,1).^3 + p(:,2).^3 + p(:,1).*(p(:,2).^2) ];
[res] = tsgEvaluate(lGrid, p);
for i = 1:4
    if (norm(res(:,i) - tv(:,i)) > 1.E-11)
        error(['Mismatch in tsgEvaluate: case 5, output ',num2str(i)]);
    end
end
tsgDeleteGrid(lGrid);

[lGrid, p] = tsgMakeGlobal('_tsgcoretests_ch', 2, 1, 'chebyshev', 'iptotal', 22);
v = [exp(-p(:,1).^2 -p(:,2).^2)];
tsgLoadValues(lGrid, v);
p = [-1.0 + 2.0 * rand(1000,2)];
v = [exp(-p(:,1).^2 -p(:,2).^2)];
[res] = tsgEvaluate(lGrid, p);
if (norm(v - res) > 1.E-9)
    error(['Mismatch in tsgEvaluate: global grid with chebyshev points']);
end
tsgDeleteGrid(lGrid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgEvaluateHierarchy()                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lGrid, p] = tsgMakeGlobal('_tsgcoretests_ml', 3, 1, 'fejer2', 'level', 4);
[V] = tsgEvaluateHierarchy(lGrid, p);
if (norm(V - eye(size(p,1))) > 1.E-11)
    error(['Mismatch in tsgEvaluateHierarchy: lagrange polynomials do not form identity']);
end

[lGrid, p] = tsgMakeSequence('_tsgcoretests_ml', 2, 1, 'leja', 'level', 2);
pnts = [0.33, 0.25; -0.27, 0.39; 0.97, -0.76; -0.44, 0.21; -0.813, 0.03; -0.666, 0.666];
tres = [ones(size(pnts, 1), 1), pnts(:,2), 0.5 * pnts(:,2) .* (pnts(:,2) - 1.0), pnts(:,1), pnts(:,1) .* pnts(:,2), 0.5 * pnts(:,1) .* (pnts(:,1) - 1.0)];
[res] = tsgEvaluateHierarchy(lGrid, pnts);
if (norm(res - tres) > 1.E-11)
    error(['Mismatch in tsgEvaluateHierarchy: sequence grid test']);
end

[lGrid, p] = tsgMakeLocalPolynomial('_tsgcoretests_ml', 2, 1, 'localp', 4, 1);
v = [exp(-p(:,1).^2 - 2.0 * p(:,2).^2)];
tsgLoadValues(lGrid, v);
pnts = [-1.0 + 2.0 * rand(13, 2)];
[tres] = tsgEvaluate(lGrid, pnts);
[mVan] = tsgEvaluateHierarchy(lGrid, pnts);
[coef] = tsgGetHCoefficients(lGrid);
res = mVan * coef;
if (norm(tres - res) > 1.E-11)
    error(['Mismatch in tsgEvaluateHierarchy: localp grid test']);
end

% this tests reading a complex matrix
[lGrid, p] = tsgMakeFourier('_tsgcoretests_ml', 2, 1, 'level', 4);
v = [exp(-p(:,1).^2 - 2.0 * p(:,2).^2)];
tsgLoadValues(lGrid, v);
pnts = [-1.0 + 2.0 * rand(13, 2)];
[tres] = tsgEvaluate(lGrid, pnts);
[mVan] = tsgEvaluateHierarchy(lGrid, pnts);
[coef] = tsgGetHCoefficients(lGrid);
res = real(mVan * coef);
if (norm(tres - res) > 1.E-11)
    error(['Mismatch in tsgEvaluateHierarchy: Fourier grid test']);
end

tsgDeleteGrid(lGrid);

[lGridA, p] = tsgMakeLocalPolynomial('_tsgcoretests_mlA', 2, 1, 'semi-localp', 5, 1);
[lGridB, p] = tsgMakeLocalPolynomial('_tsgcoretests_mlB', 2, 1, 'semi-localp', 5, 1);
v = [exp(-p(:,1).^2 - 2.0 * p(:,2).^2)];
tsgLoadValues(lGridA, v);
[mVan] = tsgEvaluateHierarchy(lGridB, p);
Coeff = mVan \ v;
tsgLoadHCoefficients(lGridB, Coeff);
pnts = [-1.0 + 2.0 * rand(13, 2)];
[tres] = tsgEvaluate(lGridA, pnts);
[res] = tsgEvaluate(lGridB, pnts);
if (norm(tres - res) > 1.E-11)
    error(['Mismatch in tsgEvaluateHierarchy: localp grid solve']);
end

[lGridA, p] = tsgMakeLocalPolynomial('_tsgcoretests_mlA', 2, 1, 'semi-localp', 5, 1, [-1 2; 7 9;]);
[lGridB, p] = tsgMakeLocalPolynomial('_tsgcoretests_mlB', 2, 1, 'semi-localp', 5, 1, [-1 2; 7 9;]);
v = [exp(-p(:,1).^2 - 2.0 * p(:,2).^2)];
tsgLoadValues(lGridA, v);
[mVan] = tsgEvaluateHierarchy(lGridB, p);
Coeff = mVan \ v;
tsgLoadHCoefficients(lGridB, Coeff);
pnts = [-1.0 + 3.0 * rand(13, 1), 7.0 + 2.0 * rand(13, 1)];
[tres] = tsgEvaluate(lGridA, pnts);
[res] = tsgEvaluate(lGridB, pnts);
if (norm(tres - res) > 1.E-11)
    error(['Mismatch in tsgEvaluateHierarchy: localp grid solve']);
end

[lGridA, p] = tsgMakeSequence('_tsgcoretests_mlA', 2, 1, 'rleja', 'level', 5, [1 2; 1 2;]);
[lGridB, p] = tsgMakeSequence('_tsgcoretests_mlB', 2, 1, 'rleja', 'level', 5, [1 2; 1 2;]);
v = [exp(-p(:,1).^2 - 2.0 * p(:,2).^2)];
tsgLoadValues(lGridA, v);
[mVan] = tsgEvaluateHierarchy(lGridB, p);
Coeff = mVan \ v;
tsgLoadHCoefficients(lGridB, Coeff);
pnts = [1.0 + rand(32, 2)];
[tres] = tsgEvaluate(lGridA, pnts);
[res] = tsgEvaluate(lGridB, pnts);
if (norm(tres - res) > 1.E-11)
    error(['Mismatch in tsgEvaluateHierarchy: sequence solve']);
end

tsgDeleteGrid(lGridA);
tsgDeleteGrid(lGridB);

% this tests the hierarchical support
[lGrid, p] = tsgMakeWavelet('_tsgcoretests_ml', 1, 1, 2, 1);
[res] = tsgGetHSupport(lGrid);
tres = [1.0, 1.0, 1.0, 1.5, 1.5, 0.75, 0.75, 0.75, 0.75]';
if (norm(tres - res) > 1.E-11)
    error(['Mismatch in tsgGetHSupport: wavelet support']);
end

tsgDeleteGrid(lGrid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgIntegrate()                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lGrid, p] = tsgMakeGlobal('_tsgcoretests_int', 1, 1, 'gauss-hermite', 'level', 2, [], [0.0, 0.0]);
v = [p.^2];
tsgLoadValues(lGrid, v)
[I] = tsgIntegrate(lGrid);
if (abs(I - pi^0.5 / 2.0) > 1.E-11)
    error('Mismatch in tsgIntegrate(): case 1');
end

[lGrid, p] = tsgMakeGlobal('_tsgcoretests_int', 1, 1, 'gauss-hermite', 'level', 2, [], [2.0, 0.0]);
v = [sqrt(2.0) * ones(size(v,1), 1)];
tsgLoadValues(lGrid, v)
[I] = tsgIntegrate(lGrid);
if (abs(I - sqrt(2.0) * pi^0.5 / 2.0) > 1.E-11)
    error('Mismatch in tsgIntegrate(): case 2');
end
tsgDeleteGrid(lGrid);

disp(['Core I/O and evaluate:    PASS']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgGetInterpolationWeights()                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lGridA, p] = tsgMakeGlobal('_tsgcoretests_IntA', 2, 1, 'fejer2', 'level', 4);
[lGridB, p] = tsgMakeGlobal('_tsgcoretests_IntB', 2, 1, 'fejer2', 'level', 4);
v = [exp(-p(:,1).^2 - 2.0 * p(:,2).^2)];
tsgLoadValues(lGridA, v);
pnts = [-1.0 + 2.0 * rand(32, 2)];
[tres] = tsgEvaluate(lGridA, pnts);
[A] = tsgGetInterpolationWeights(lGridB, pnts);
res = A * v;
if (norm(tres - res) > 1.E-11)
    error(['Mismatch in tsgGetInterpolationWeights: global case']);
end

[lGridA, p] = tsgMakeSequence('_tsgcoretests_IntA', 2, 1, 'min-delta', 'level', 7);
[lGridB, p] = tsgMakeSequence('_tsgcoretests_IntB', 2, 1, 'min-delta', 'level', 7);
v = [exp(-p(:,1).^2 - 2.0 * p(:,2).^2)];
tsgLoadValues(lGridA, v);
pnts = [-1.0 + 2.0 * rand(32, 2)];
[tres] = tsgEvaluate(lGridA, pnts);
[A] = tsgGetInterpolationWeights(lGridB, pnts);
res = A * v;
if (norm(tres - res) > 1.E-11)
    error(['Mismatch in tsgGetInterpolationWeights: global case']);
end

tsgDeleteGrid(lGridA);
tsgDeleteGrid(lGridB);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgEstimateAnisotropicCoefficients()         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lGrid, p] = tsgMakeGlobal('_tsgcoretests_ans', 2, 1, 'rleja', 'level', 9);
v = [exp(p(:,1) + p(:,2).^2)];
tsgLoadValues(lGrid, v);
[c] = tsgEstimateAnisotropicCoefficients(lGrid, 'iptotal');
if (abs(c(1) / c(2) - 2.0) > 0.2)
    error('Mismatch in tsgEstimateAnisotropicCoefficients(): total degree');
end
[c] = tsgEstimateAnisotropicCoefficients(lGrid, 'ipcurved');
if (length(c) ~= 4)
    error('Mismatch in tsgEstimateAnisotropicCoefficients(): curved dimensions');
end
if ((abs(c(1) / c(2) - 2.0) > 0.2) || (c(3) > 0.0) || (c(4) > 0.0))
    error('Mismatch in tsgEstimateAnisotropicCoefficients(): curved');
end
tsgDeleteGrid(lGrid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgRefineAnisotropic()                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lGrid, p] = tsgMakeSequence('_tsgcoretests_refA', 3, 1, 'leja', 'level', 3, [], [], [], [], [3, 2, 1]);
if ((sum((abs(p(:,1) - sqrt(1.0/3.0)) < 0.0001)) == 0) || (sum((abs(p(:,2) - sqrt(1.0/3.0)) < 0.0001)) > 0) || (sum((abs(p(:,3) - sqrt(1.0/3.0)) < 0.0001)) > 0))
    error('Mismatch in tsgRefineAnisotropic(): limits in make');
end
v = [exp(-p(:,1).^2 - p(:,2).^2)];
tsgLoadValues(lGrid, v);
tsgRefineAnisotropic(lGrid, 'iptotal', 5, 0);
[p] = tsgGetNeededPoints(lGrid);
if (size(p, 1) == 0)
    error('Mismatch in tsgRefineAnisotropic(): did not refine');
end
if ((sum((abs(p(:,2) - sqrt(1.0/3.0)) < 0.0001)) > 0) || (sum((abs(p(:,3) - sqrt(1.0/3.0)) < 0.0001)) > 0))
    error('Mismatch in tsgRefineAnisotropic(): limits refine using existing limits');
end
tsgRefineAnisotropic(lGrid, 'iptotal', 10, 0, [3, 2, 2]);
[p] = tsgGetNeededPoints(lGrid);
if (size(p, 1) == 0)
    error('Mismatch in tsgRefineAnisotropic(): did not refine');
end
if ((sum((abs(p(:,2) - sqrt(1.0/3.0)) < 0.0001)) > 0) || (sum((abs(p(:,3) - 1.0) < 0.0001)) == 0) || (sum((abs(p(:,3) - sqrt(1.0/3.0)) < 0.0001)) > 0))
    error('Mismatch in tsgRefineAnisotropic(): limits refine using new limits');
end
tsgDeleteGrid(lGrid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgRefineSurplus()                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lGrid, p] = tsgMakeLocalPolynomial('_tsgcoretests_refS', 3, 1, 'localp', 3, 1, [], [], [], [1, 2, 3]);
if ((sum((abs(p(:,1) - 0.5) < 0.0001)) > 0) || (sum((abs(p(:,2) - 0.25) < 0.0001)) > 0))
    error('Mismatch in tsgRefineSurplus(): limits in make');
end
v = [exp(-p(:,1).^2 - p(:,2).^2)];
tsgLoadValues(lGrid, v);
tsgRefineSurplus(lGrid, 1.E-8, 'classic', 0);
[p] = tsgGetNeededPoints(lGrid);
if (size(p, 1) == 0)
    error('Mismatch in tsgRefineSurplus(): did not refine local polynomial');
end
if ((sum((abs(p(:,1) - 0.5) < 0.0001)) > 0) || (sum((abs(p(:,2) - 0.25) < 0.0001)) > 0))
    error('Mismatch in tsgRefineSurplus(): limits refine using existing limits');
end
tsgRefineSurplus(lGrid, 1.E-8, 'classic', 0, [2, 2, 3]);
[p] = tsgGetNeededPoints(lGrid);
if (size(p, 1) == 0)
    error('Mismatch in tsgRefineSurplus(): did not refine on second pass');
end
if ((sum((abs(p(:,1) - 0.5) < 0.0001)) == 0) || (sum((abs(p(:,1) - 0.25) < 0.0001)) > 0) || (sum((abs(p(:,2) - 0.25) < 0.0001)) > 0))
    error('Mismatch in tsgRefineSurplus(): limits refine using new limits');
end
tsgDeleteGrid(lGrid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgCancelRefine()                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lGrid, p] = tsgMakeLocalPolynomial('_tsgcoretests_cref', 3, 1, 'localp', 4, 2);
v = [exp(-p(:,1).^2 -0.5 * p(:,2).^2 -2.0 * p(:,3).^2)];
tsgLoadValues(lGrid, v);
[p] = tsgGetNeededPoints(lGrid);
if (max(size(p)) > 0)
    error('Mismatch in cancel refine: did not load values');
end
tsgRefineSurplus(lGrid, 1.E-4, 'direction');
[p] = tsgGetNeededPoints(lGrid);
if (max(size(p)) == 0)
    error('Mismatch in cancel refine: did not set refinement at output -1');
end
tsgCancelRefine(lGrid);
[p] = tsgGetNeededPoints(lGrid);
if (max(size(p)) > 0)
    error('Mismatch in cancel refine: did not cancel the refinement');
end
tsgDeleteGrid(lGrid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgMergeRefine()                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lGridA, p] = tsgMakeGlobal('_tsgcoretests_refA', 2, 1, 'fejer2', 'level', 4);
[lGridB, p] = tsgMakeGlobal('_tsgcoretests_refB', 2, 1, 'fejer2', 'level', 4);
v = [exp(-p(:,1).^2 -p(:,2).^2)];
tsgLoadValues(lGridA, v);
tsgLoadValues(lGridB, v);
tsgRefineAnisotropic(lGridA, 'iptotal', 10, 0);
tsgRefineAnisotropic(lGridB, 'iptotal', 10, 0);
[p] = tsgGetNeededPoints(lGridA);
v = [exp(-p(:,1).^2 -p(:,2).^2)];
tsgLoadValues(lGridA, v);
tsgMergeRefine(lGridB)
[p1] = tsgGetPoints(lGridA);
[p2] = tsgGetPoints(lGridB);
if (norm(p1 - p2) > 1.E-11)
    error('Mismatch in tsgMergeRefine(): case 2, tsgGetPoints()')
end
p = [-1.0 + 2.0 * rand(20, 2)];
[vB] = tsgEvaluate(lGridB, p);
if (norm(vB) > 1.E-11)
    error('Mismatch in tsgMergeRefine(): case 3, tsgEvaluate() not zero');
end
[p] = tsgGetPoints(lGridB);
v = [exp(-p(:,1).^2 -p(:,2).^2)];
tsgLoadValues(lGridB, v);
p = [-1.0 + 2.0 * rand(30, 2)];
[vA] = tsgEvaluate(lGridA, p);
[vB] = tsgEvaluate(lGridB, p);
if (norm(vA - vB) > 1.E-11)
    error('Mismatch in tsgMergeRefine(): case 3, tsgEvaluate() not equal');
end
tsgDeleteGrid(lGridA);
tsgDeleteGrid(lGridB);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgGetHCoefficients()                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% covered in tsgEvaluateHierarchy()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgLoadHCoefficients()                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% covered in tsgEvaluateHierarchy()

disp(['Refinement functions:     PASS']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     tsgReloadGrid()                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lGrid, p] = tsgMakeGlobal('_tsgcoretests_ch', 3, 7, 'chebyshev', 'iptotal', 5);
[lGrid2] = tsgReloadGrid('_tsgcoretests_ch');
if (lGrid2.sName ~= '_tsgcoretests_ch')
    error('Mismatch in tsgReloadGrid() could not reload grid: sName');
end
if (lGrid2.sType ~= 'Global')
    error('Mismatch in tsgReloadGrid() could not reload grid: sType');
end
if ((lGrid2.iDim ~= 3) || (lGrid2.iOut ~= 7))
    error('Mismatch in tsgReloadGrid() could not reload grid: iDim and iOut');
end
tsgDeleteGrid(lGrid);

disp(['Utility functions:        PASS']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['']);
disp(['All Tasmanian Tests Completed Successfully']);

end
