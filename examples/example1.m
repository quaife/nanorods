%close all

prams.N = 64; % points per body
prams.nv = 20; % number of bodies
prams.T = 10/100; % time horizon
prams.m = 100/100; % number of time steps
prams.lengths = 4.5*ones(1, prams.nv);
prams.widths = 1*ones(1,prams.nv);
prams.order = 4;

options.farField = 'poiseuille';
options.usePlot = true;
options.axis = [-20 20 -5 5];
options.saveData = true;
options.dataFile = 'rectangular_fibers_poiseuille';
options.append = false;
options.inear = true;
options.usePreco = true;
options.ifmm = true;

[options,prams] = initRigid2D(options,prams);

%% staggerd grid
rown = 5;
coln = prams.nv/(2*rown);

x = linspace(0, 3*max(prams.widths)*rown, rown);
y = linspace(0, 0.75*max(prams.lengths)*coln, coln);


[X1, Y1] = meshgrid(x,y);
[X2, Y2] = meshgrid(x - 2*max(prams.widths), y - max(prams.lengths)/2);

coeffr = 0;
xc = [[X1(:)', X2(:)'] + coeffr*(1 - 2*rand(1,prams.nv)); ...
    [Y1(:)' - coln/2*max(prams.lengths), Y2(:)'- coln/2*max(prams.lengths)] + coeffr*(1 - 2*rand(1,prams.nv))];
tau = pi/2*ones(1,prams.nv) + 2*coeffr*(1-2*rand(1,prams.nv));

Xfinal = rigid2D(options, prams, xc, tau);

%pp = post([options.dataFile,'.dat']);
%pp.animated_gif('rectangular_fibers_poiseuille.gif', 1, [])
%stats = pp.calculate_stats(1:prams.nv);
