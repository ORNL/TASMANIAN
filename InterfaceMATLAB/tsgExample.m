function tsgExample(bFast)
%
% tsgExample()
%
% this is example source code on how to call the different functions
% this does the exact same thing as the C++ example
%

%  bFast: if this exists, then only the first two tests are executed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EXAMPLE 1: 
%
% integrate: f(x,y) = exp(-x^2) * cos(y) over [-1,1] x [-1,1]
% using classical Smolyak grid with Clenshaw-Curtis points and weights
%

dim = 2;
level = 6;
order = 0; % not used by Clenshaw-Curtis rule
[weights, points] = tsgMakeQuadrature(dim, 'clenshaw-curtis', 'level', level, order);

I = weights' * (exp(-points(:,1).^2) .* cos(points(:,2))); % I is the approximated quadrature
E = 2.513723354063905e+00; % E is the "exact" solution computed to 16 decimal places
E = abs(I - E);
%E = abs(I - quad(@(x)(exp(-x.^2)), -1, 1, 1.E-14)*quad(@(x)(cos(x)), -1, 1, 1.E-14));

disp(['----------------------------------------------------------------------------']);
disp([' Example 1:  integrate f(x,y) = exp(-x^2) * cos(y), using clenshaw-curtis level nodes']);
disp(['    at level ',num2str(level)]);
disp(['    grid has ',num2str(size(points,1)),' points']);
disp(['    integral: ',num2str(I,16)]);
disp(['    error: ',num2str(E,16)]);
disp([' ']);

level = 7;
[weights, points] = tsgMakeQuadrature(dim, 'clenshaw-curtis', 'level', level, order);

I = weights' * (exp(-points(:,1).^2) .* cos(points(:,2))); % I is the approximated quadrature
E = 2.513723354063905e+00; % E is the "exact" solution computed to 16 decimal places
E = abs(I - E);

disp(['    at level ',num2str(level)]);
disp(['    grid has ',num2str(size(points,1)),' points']);
disp(['    integral: ',num2str(I,16)]);
disp(['    error: ',num2str(E,16),' (rounded to 14 decimal places)']);
disp([' ']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EXAMPLE 2: 
%
% integrate: f(x,y) = exp(-x^2) * cos(y) over (x,y) in [-5,5] x [-2,3]
% using Gauss-Patterson rules chosen to integrate exactly polynomials of
% total degree up to degree specified by prec
%

dim = 2;
prec = 20;
domain = [-5, 5; -2, 3];
order = 0; % not used by Gauss-Patterson rule
[weights, points] = tsgMakeQuadrature(dim, 'gauss-patterson', 'qptotal', prec, order, domain);

I = weights' * (exp(-points(:,1).^2) .* cos(points(:,2))); % I is the approximated quadrature
E = 1.861816427518323e+00; % E is the "exact" solution computed to 16 decimal places
E = abs(I - E);
%E = abs(I - quad(@(x)(exp(-x.^2)), -1, 1, 1.E-14)*quad(@(x)(cos(x)), -1, 1, 1.E-14));

disp(['----------------------------------------------------------------------------']);
disp([' Example 2: integrate f(x,y) = exp(-x^2) * cos(y) over [-5,5] x [-2,3] using  Gauss-Patterson nodes']);
disp(['    at precision ',num2str(prec)]);
disp(['    grid has ',num2str(size(points,1)),' points']);
disp(['    integral: ',num2str(I,16)]);
disp(['    error: ',num2str(E,16)]);
disp([' ']);

prec = 40;
[weights, points] = tsgMakeQuadrature(dim, 'gauss-patterson', 'qptotal', prec, order, domain);

I = weights' * (exp(-points(:,1).^2) .* cos(points(:,2))); % I is the approximated quadrature
E = 1.861816427518323e+00; % E is the "exact" solution computed to 16 decimal places
E = abs(I - E);

disp(['    at precision ',num2str(prec)]);
disp(['    grid has ',num2str(size(points,1)),' points']);
disp(['    integral: ',num2str(I,16)]);
disp(['    error: ',num2str(E,16)]);
disp([' ']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (exist('bFast')) % this is used for automated testing, ignore
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EXAMPLE 3: 
%
% integrate: f(x,y) = exp(-x^2) * cos(y) over (x,y) in [-5,5] x [-2,3]
% using Gauss-Patterson, Clenshaw-Curtis and Gauss-Legendre rules and
% compare the results
%

dim = 2;
domain = [-5, 5; -2, 3];
order = 0;
E = 1.861816427518323e+00; % E is the "exact" solution computed to 16 decimal places

disp(['----------------------------------------------------------------------------']);
disp([' Example 3: integrate f(x,y) = exp(-x^2) * cos(y) over [-5,5] x [-2,3] using different rules']);
disp([' ']);
disp(['             Clenshaw-Curtis         Gauss-Legendre          Gauss-Patterson']);
disp([' precision   nodes       error       nodes       error       nodes       error']);

for prec = 9:4:30
    tt = num2str(prec);
    ss = [blanks(6 - length(tt)),tt];
    
    [weights, points] = tsgMakeQuadrature(dim, 'clenshaw-curtis', 'qptotal', prec, order, domain);
    I = weights' * (exp(-points(:,1).^2) .* cos(points(:,2)));
    tt = num2str(size(points,1));
    ss = [ss,'      ',blanks(5 - length(tt)),tt];
    tt = num2str(abs(I-E),5);
    ss = [ss,' ',blanks(12 - length(tt)),tt];
    
    ss = [ss,blanks(30 - length(ss))];
    
    [weights, points] = tsgMakeQuadrature(dim, 'gauss-legendre', 'qptotal', prec, order, domain);
    I = weights' * (exp(-points(:,1).^2) .* cos(points(:,2)));
    tt = num2str(size(points,1));
    ss = [ss,'      ',blanks(5 - length(tt)),tt];
    tt = num2str(abs(I-E),5);
    ss = [ss,' ',blanks(12 - length(tt)),tt];
    
    ss = [ss,blanks(60 - length(ss))];
    
    [weights, points] = tsgMakeQuadrature(dim, 'gauss-patterson', 'qptotal', prec, order, domain);
    I = weights' * (exp(-points(:,1).^2) .* cos(points(:,2)));
    tt = num2str(size(points,1));
    ss = [ss,blanks(5 - length(tt)),tt];
    tt = num2str(abs(I-E),5);
    ss = [ss,' ',blanks(12 - length(tt)),tt];
    
    disp(ss);
end
disp([' ']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EXAMPLE 4:
%
% interpolate: f(x,y) = exp(-x^2) * cos(y)
% over [2,3] x [2,3]
% with a rule that exactly interpolates polynomials of total degree up to 
% degree specified by prec
%
% NOTE: any grid with name '_tsgExample4' will be overwritten by the
% tsgMakeGlobal() command
%

dim = 2;
outs = 1;
prec = 6;
[lGrid, points] = tsgMakeGlobal('_tsgExample4', dim, outs, 'clenshaw-curtis', 'iptotal', prec, [2 3; 2 3]);

vals = (exp(-points(:,1).^2) .* cos(points(:,2)));
tsgLoadValues(lGrid, vals);

[res] = tsgEvaluate(lGrid, [2.3, 2.7]);

disp(['----------------------------------------------------------------------------']);
disp([' Example 4: interpolate f(x,y) = exp(-x^2) * cos(y) over [2,3] x [2,3]']);
disp(['    using clenshaw-curtis iptotal rule']);
disp(['    using polynomials of total degree up to ',num2str(prec)]);
disp(['    grid has ',num2str(size(points,1)),' points']);
disp(['    interpolant at (2.3,2.7): ',num2str(res,16)]);
disp(['    error: ',num2str(abs(res - exp(-2.3^2) * cos(2.7)),16)]);
disp([' ']);

prec = 12;
[lGrid, points] = tsgMakeGlobal('_tsgExample4', dim, outs, 'clenshaw-curtis', 'iptotal', prec, [2 3; 2 3]);

vals = (exp(-points(:,1).^2) .* cos(points(:,2)));
tsgLoadValues(lGrid, vals);

[res] = tsgEvaluate(lGrid, [2.3, 2.7]);

disp(['    using polynomials of total degree up to ',num2str(prec)]);
disp(['    grid has ',num2str(size(points,1)),' points']);
disp(['    interpolant at (2.3,2.7): ',num2str(res,16)]);
disp(['    error: ',num2str(abs(res - exp(-2.3^2) * cos(2.7)),16)]);
disp([' ']);
tsgDeleteGrid(lGrid); % clear all temporary used files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EXAMPLE 5:
%
% interpolate: f(x1,x2,x3,x4) = exp(-x1^2) * cos(x2) * exp(-x3^2) * cos(x4)
% with Global and Sequence Leja rules
%

dim = 4;
outs = 1;
prec = 15;

tic;
[lGrid, points] = tsgMakeGlobal('_tsgExample5', dim, outs, 'leja', 'iptotal', prec);
gstage1 = toc;

disp(['----------------------------------------------------------------------------']);
disp([' Example 5: interpolate f(x1,x2,x3,x4) = exp(-x1^2) * cos(x2) * exp(-x3^2) * cos(x4)']);
disp(['    comparign the performance of Global and Sequence grids with leja nodes']);
disp(['    using polynomials of total degree up to ',num2str(prec)]);
disp(['    grid has ',num2str(size(points,1)),' points']);
disp(['    both grids are evaluated at 1000 random points']);
disp([' ']);


vals = (exp(-points(:,1).^2) .* cos(points(:,2))) .* (exp(-points(:,3).^2) .* cos(points(:,4)));

tic;
tsgLoadValues(lGrid, vals);
gstage2 = toc;

pnts = [-1 + 2 * rand(1000, 4)];
tres = (exp(-pnts(:,1).^2) .* cos(pnts(:,2))) .* (exp(-pnts(:,3).^2) .* cos(pnts(:,4)));

tic;
[res] = tsgEvaluate(lGrid, pnts);
gstage3 = toc;

gerr = max(abs(res - tres));


tic;
[lGrid, points] = tsgMakeSequence('_tsgExample5', dim, outs, 'leja', 'iptotal', prec);
sstage1 = toc;

vals = (exp(-points(:,1).^2) .* cos(points(:,2))) .* (exp(-points(:,3).^2) .* cos(points(:,4)));

tic;
tsgLoadValues(lGrid, vals);
sstage2 = toc;

tic;
[sres] = tsgEvaluate(lGrid, pnts);
sstage3 = toc;

serr = max(abs(sres - tres));


disp(['    Stage         Global Grid       Sequence Grid']);
disp(['    make grid       ',num2str(gstage1),'          ',num2str(sstage1),'    seconds']);
disp(['    load values     ',num2str(gstage2),'          ',num2str(sstage2),'    seconds']);
disp(['    evaluate        ',num2str(gstage3),'          ',num2str(sstage3),'    seconds']);
disp(['    error           ',num2str(gerr),'       ',num2str(serr)]);
disp(['    The difference between the two approximations is: ',num2str(max(abs(res - sres)))]);
disp(['  NOTE: the MATLAB interface has the additional overhead of reading/writing files to the work folder']);
disp(['        the results here are not the best comparison between sequence and global rules']);
disp([' ']);
tsgDeleteGrid(lGrid); % clear all temporary used files

%
% I tried to make an example here with more outputs, but the cost of
% tsgLoadValues() is dominated by the reading/writing of the files and is
% hence not a good comparison. See the C++ example.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EXAMPLE 6:
%
% interpolate: f(x,y) = exp(-x^2) * cos(y) 
% using different refinement schemes
%

dim = 2;
outs = 1;
prec = 3;

[lGrid1, points] = tsgMakeGlobal('_tsgExample6a', dim, outs, 'leja', 'iptotal', prec);
[lGrid2, points] = tsgMakeGlobal('_tsgExample6b', dim, outs, 'leja', 'iptotal', prec);
[lGrid3, points] = tsgMakeGlobal('_tsgExample6c', dim, outs, 'leja', 'iptotal', prec);

vals = (exp(-points(:,1).^2) .* cos(points(:,2)));

tsgLoadValues(lGrid1, vals);
tsgLoadValues(lGrid2, vals);
tsgLoadValues(lGrid3, vals);

pnts = [-1 + 2 * rand(1000, 2)];
tres = exp(-pnts(:,1).^2) .* cos(pnts(:,2));

disp(['----------------------------------------------------------------------------']);
disp([' Example 6: interpolate: f(x,y) = exp(-x^2) * cos(y) ']);
disp(['    using leja nodes and different refinement schemes']);
disp(['    the error is estimated as the maximum from 1000 random points']);
disp([' ']); % clear all temporary used files

disp(['             Total Degree            Curved                  Surplus']);
disp([' iteration   nodes       error       nodes       error       nodes       error']);

nump1 = size(points, 1);
nump2 = size(points, 1);
nump3 = size(points, 1);

for iI = 1:10
    tt = num2str(iI);
    ss = [blanks(6 - length(tt)),tt];
    
    [points] = tsgRefineAnisotropic(lGrid1, 'iptotal', 10);
    vals = (exp(-points(:,1).^2) .* cos(points(:,2)));
    tsgLoadValues(lGrid1, vals);
    nump1 = nump1 + size(points, 1);
    tt = num2str(nump1);
    ss = [ss,'      ',blanks(5 - length(tt)),tt];
    [res] = tsgEvaluate(lGrid1, pnts);
    tt = num2str(max(abs(res - tres)),5);
    ss = [ss,' ',blanks(12 - length(tt)),tt];
    
    ss = [ss,blanks(30 - length(ss))];
    
    [points] = tsgRefineAnisotropic(lGrid2, 'ipcurved', 10);
    vals = (exp(-points(:,1).^2) .* cos(points(:,2)));
    tsgLoadValues(lGrid2, vals);
    nump2 = nump2 + size(points, 1);
    tt = num2str(nump2);
    ss = [ss,'      ',blanks(5 - length(tt)),tt];
    [res] = tsgEvaluate(lGrid2, pnts);
    tt = num2str(max(abs(res - tres)),5);
    ss = [ss,' ',blanks(12 - length(tt)),tt];
    
    ss = [ss,blanks(60 - length(ss))];
    
    [points] = tsgRefineSurplus(lGrid3, 1.E-10);
    vals = (exp(-points(:,1).^2) .* cos(points(:,2)));
    tsgLoadValues(lGrid3, vals);
    nump3 = nump3 + size(points, 1);
    tt = num2str(nump3);
    ss = [ss,blanks(5 - length(tt)),tt];
    [res] = tsgEvaluate(lGrid3, pnts);
    tt = num2str(max(abs(res - tres)),5);
    ss = [ss,' ',blanks(12 - length(tt)),tt];
    
    disp(ss);
end
disp([' ']);

tsgDeleteGrid(lGrid1);
tsgDeleteGrid(lGrid2);
tsgDeleteGrid(lGrid3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EXAMPLE 7:
%
% interpolate: f(x,y) = exp(-x^2) * cos(y) 
% using different local polynomial rules
%

dim = 2;
outs = 1;
prec = 7;

disp(['----------------------------------------------------------------------------']);
disp([' Example 7: interpolate: f(x,y) = exp(-x^2) * cos(y) ']);
disp(['    using localp and semi-localp rules with depth ',num2str(prec)]);
disp(['    the error is estimated as the maximum from 1000 random points']);
disp([' ']); % clear all temporary used files

[lGrid, points] = tsgMakeLocalPolynomial('_tsgExample7', dim, outs, 'localp', prec, 2);
vals = (exp(-points(:,1).^2) .* cos(points(:,2)));
tsgLoadValues(lGrid, vals);
pnts = [-1 + 2 * rand(1000, 2)];
tres = exp(-pnts(:,1).^2) .* cos(pnts(:,2));
[res] = tsgEvaluate(lGrid, pnts);

disp(['   Number of points: ',num2str(size(points, 1))]);
disp(['   Error for      localp: ',num2str(max(abs(res - tres)))]);

[lGrid, points] = tsgMakeLocalPolynomial('_tsgExample7', dim, outs, 'semi-localp', prec, 2);
vals = (exp(-points(:,1).^2) .* cos(points(:,2)));
tsgLoadValues(lGrid, vals);
[res] = tsgEvaluate(lGrid, pnts);

disp(['   Error for semi-localp: ',num2str(max(abs(res - tres)))]);
disp([' Note: semi-localp wins this competition because the function is very smooth']);
disp([' ']);

tsgDeleteGrid(lGrid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EXAMPLE 8:
%
% interpolate f(x,y) = cos(0.5 * pi * x) * cos(0.5 * pi * y)
% using different local polynomial rules
%

dim = 2;
outs = 1;
prec = 7;

disp(['----------------------------------------------------------------------------']);
disp([' Example 8: interpolate f(x,y) = cos(0.5 * pi * x) * cos(0.5 * pi * y) ']);
disp(['    using localp and localp-zero rules with depth ',num2str(prec)]);
disp(['    the error is estimated as the maximum from 1000 random points']);
disp([' ']); 

[lGrid, points] = tsgMakeLocalPolynomial('_tsgExample8', dim, outs, 'localp', prec, 2);
vals = (cos(0.5 * pi * points(:,1)) .* cos(0.5 * pi * points(:,2)));
tsgLoadValues(lGrid, vals);
pnts = [-1 + 2 * rand(1000, 2)];
tres = (cos(0.5 * pi * pnts(:,1)) .* cos(0.5 * pi * pnts(:,2)));
[res] = tsgEvaluate(lGrid, pnts);

disp(['   localp       Number of points: ',num2str(size(points, 1)),'  Error: ',num2str(max(abs(res - tres)))]);

[lGrid, points] = tsgMakeLocalPolynomial('_tsgExample8', dim, outs, 'localp-zero', prec-1, 2);
vals = (cos(0.5 * pi * points(:,1)) .* cos(0.5 * pi * points(:,2)));
tsgLoadValues(lGrid, vals);
[res] = tsgEvaluate(lGrid, pnts);

disp(['   localp-zero  Number of points: ',num2str(size(points, 1)),'  Error: ',num2str(max(abs(res - tres)))]);
disp([' Note: localp-zero wins this competition because the function is zero at the boundary']);
disp([' ']);

tsgDeleteGrid(lGrid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EXAMPLE 9:
%
% interpolate f(x,y) = exp(- x) / (1 + 100 * exp(- 10 * y))
% using different local refinement schemes
%

dim = 2;
outs = 1;
prec = 2;
tol = 1.E-5;

disp(['----------------------------------------------------------------------------']);
disp([' Example 9: interpolate f(x,y) = exp(- x) / (1 + 100 * exp(- 10 * y)) ']);
disp(['    the error is estimated as the maximum from 1000 random points']);
disp(['    tolerance is set at 1.E-5 and maximal order polynomials are used']);
disp([' ']);

[lGrid1, points] = tsgMakeLocalPolynomial('_tsgExample9a', dim, outs, 'localp', prec, -1);
vals = exp(-points(:,1)) ./ (1 + 100 * exp(-10 * points(:,2)));
tsgLoadValues(lGrid1, vals);
pnts = [-1 + 2 * rand(1000, 2)];
tres = exp(-pnts(:,1)) ./ (1 + 100 * exp(-10 * pnts(:,2)));

[lGrid2, points] = tsgMakeLocalPolynomial('_tsgExample9b', dim, outs, 'localp', prec, -1);
tsgLoadValues(lGrid2, vals);

disp(['             Classic                 FDS']);
disp([' iteration   nodes       error       nodes       error']);

nump1 = size(points, 1);
nump2 = size(points, 1);

for iI = 1:7
    tt = num2str(iI);
    ss = [blanks(6 - length(tt)),tt];
    
    [points] = tsgRefineSurplus(lGrid1, tol, 'classic');
    vals = exp(-points(:,1)) ./ (1 + 100 * exp(-10 * points(:,2)));
    tsgLoadValues(lGrid1, vals);
    nump1 = nump1 + size(points, 1);
    tt = num2str(nump1);
    ss = [ss,'      ',blanks(5 - length(tt)),tt];
    [res] = tsgEvaluate(lGrid1, pnts);
    tt = num2str(max(abs(res - tres)),5);
    ss = [ss,' ',blanks(12 - length(tt)),tt];
    
    ss = [ss,blanks(30 - length(ss))];
    
    [points] = tsgRefineSurplus(lGrid2, tol, 'fds');
    vals = exp(-points(:,1)) ./ (1 + 100 * exp(-10 * points(:,2)));
    tsgLoadValues(lGrid2, vals);
    nump2 = nump2 + size(points, 1);
    tt = num2str(nump2);
    ss = [ss,'      ',blanks(5 - length(tt)),tt];
    [res] = tsgEvaluate(lGrid2, pnts);
    tt = num2str(max(abs(res - tres)),5);
    ss = [ss,' ',blanks(12 - length(tt)),tt];
    
    disp(ss);
end
disp([' ']);

tsgDeleteGrid(lGrid1);
tsgDeleteGrid(lGrid2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EXAMPLE 10:
%
% interpolate f(x,y) = exp(- x) / (1 + 100 * exp(- 10 * y))
% using local polynomails and wavelets
%

dim = 2;
outs = 1;
prec = 1;
tol = 1.E-5;

disp(['----------------------------------------------------------------------------']);
disp([' Example 10: interpolate f(x,y) = exp(- x) / (1 + 100 * exp(- 10 * y)) ']);
disp(['    using local polynomials and wavelets']);
disp(['    the error is estimated as the maximum from 1000 random points']);
disp([' ']);

[lGrid1, points] = tsgMakeLocalPolynomial('_tsgExample10a', dim, outs, 'localp', prec+2, 1);
vals = exp(-points(:,1)) ./ (1 + 100 * exp(-10 * points(:,2)));
tsgLoadValues(lGrid1, vals);
pnts = [-1 + 2 * rand(1000, 2)];
tres = exp(-pnts(:,1)) ./ (1 + 100 * exp(-10 * pnts(:,2)));

[lGrid2, points] = tsgMakeWavelet('_tsgExample10b', dim, outs, prec, 1);
vals = exp(-points(:,1)) ./ (1 + 100 * exp(-10 * points(:,2)));
tsgLoadValues(lGrid2, vals);

disp(['             Polynomials             Wavelets']);
disp([' iteration   nodes       error       nodes       error']);

nump1 = size(points, 1);
nump2 = size(points, 1);

for iI = 1:8
    tt = num2str(iI);
    ss = [blanks(6 - length(tt)),tt];
    
    [points] = tsgRefineSurplus(lGrid1, tol, 'fds');
    vals = exp(-points(:,1)) ./ (1 + 100 * exp(-10 * points(:,2)));
    tsgLoadValues(lGrid1, vals);
    nump1 = nump1 + size(points, 1);
    tt = num2str(nump1);
    ss = [ss,'      ',blanks(5 - length(tt)),tt];
    [res] = tsgEvaluate(lGrid1, pnts);
    tt = num2str(max(abs(res - tres)),5);
    ss = [ss,' ',blanks(12 - length(tt)),tt];
    
    ss = [ss,blanks(30 - length(ss))];
    
    [points] = tsgRefineSurplus(lGrid2, tol, 'fds');
    vals = exp(-points(:,1)) ./ (1 + 100 * exp(-10 * points(:,2)));
    tsgLoadValues(lGrid2, vals);
    nump2 = nump2 + size(points, 1);
    tt = num2str(nump2);
    ss = [ss,'      ',blanks(5 - length(tt)),tt];
    [res] = tsgEvaluate(lGrid2, pnts);
    tt = num2str(max(abs(res - tres)),5);
    ss = [ss,' ',blanks(12 - length(tt)),tt];
    
    disp(ss);
end
disp([' Note: wavelets have a larger Lebesgue constant and thus wavelets are not always better than polynomials.']);
disp([' ']);

tsgDeleteGrid(lGrid1);
tsgDeleteGrid(lGrid2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EXAMPLE 11: 
%
% interpolate: f(x,y,z) = 1/((1+4x^2)*(1+5y^2)*(1+6z^2))
% using classical and conformal transformation
%
dim = 3;
prec = 12;
outs = 1;
[lGrid, points] = tsgMakeGlobal('_tsgExample11a', dim, outs, 'clenshaw-curtis', 'iptotal', prec, [], [], [], []);
vals = 1.0 ./ ((1.0 + 4.0*points(:,1).^2).*(1.0 + 5.0*points(:,2).^2).*(1.0 + 6.0*points(:,3).^2));
tsgLoadValues(lGrid, vals);
pnts = [-1 + 2 * rand(1000, 3)];
tres = 1.0 ./ ((1.0 + 4.0*pnts(:,1).^2).*(1.0 + 5.0*pnts(:,2).^2).*(1.0 + 6.0*pnts(:,3).^2));
[res] = tsgEvaluate(lGrid, pnts);
err1 = max(abs(tres - res));
nump1 = size(points, 1);

[lGrid, points] = tsgMakeGlobal('_tsgExample11a', dim, outs, 'clenshaw-curtis', 'iptotal', prec, [], [], [], [], 'asin', [4, 4, 4]);
vals = 1.0 ./ ((1.0 + 4.0*points(:,1).^2).*(1.0 + 5.0*points(:,2).^2).*(1.0 + 6.0*points(:,3).^2));
tsgLoadValues(lGrid, vals);
[res] = tsgEvaluate(lGrid, pnts);
err2 = max(abs(tres - res));

[lGrid, points] = tsgMakeLocalPolynomial('_tsgExample11a', dim, outs, 'localp', prec-4, 2);
vals = 1.0 ./ ((1.0 + 4.0*points(:,1).^2).*(1.0 + 5.0*points(:,2).^2).*(1.0 + 6.0*points(:,3).^2));
tsgLoadValues(lGrid, vals);
err3 = max(abs(tres - res));

[lGrid, points] = tsgMakeLocalPolynomial('_tsgExample11a', dim, outs, 'localp', prec-4, 2, [], 'asin', [4, 4, 4]);
vals = 1.0 ./ ((1.0 + 4.0*points(:,1).^2).*(1.0 + 5.0*points(:,2).^2).*(1.0 + 6.0*points(:,3).^2));
tsgLoadValues(lGrid, vals);
[res] = tsgEvaluate(lGrid, pnts);
err4 = max(abs(tres - res));
nump2 = size(points, 1);

disp(['----------------------------------------------------------------------------']);
disp([' Example 11: interpolate f(x,y,z) = 1/((1+4x^2)*(1+5y^2)*(1+6z^2))']);
disp(['             using conformal transformation']);
disp(['             the error is estimated as the maximum from 1000 random points']);
disp(['']);

disp([' Grid Type    nodes     error regular   error conformal']);
disp([' Global         ', num2str(nump1,'%5d'),'         ', num2str(err1, '%1.3e'), '        ', num2str(err2, '%1.3e')]);
disp([' Localp        ',  num2str(nump2,'%5d'),'         ', num2str(err3, '%1.3e'), '        ', num2str(err4, '%1.3e')]);

disp([' Note: conformal maps address specific problems with the region of analyticity of a function']);
disp(['       the map can accelerate or slow down convergence depending on the problem']);
disp([' ']);

tsgDeleteGrid(lGrid);

end

