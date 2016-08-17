close all

prams.N = 32; % points per body
prams.Nbd = 192; %points on solid wall
prams.nv = 1; % number of bodies
prams.nbd = 1; %number of walls
prams.T = 1; % time horizon
prams.m = 10; % number of time steps
prams.lengths = 4.5*ones(1, prams.nv);
prams.widths = 1*ones(1,prams.nv);
prams.order = 4;

options.farField = 'couette';
options.saveData = true;
options.fileBase = 'confined_flow';
options.append = false;
options.inear = true;
options.usePreco = false;
options.ifmm = true;
options.verbose = true;
options.profile = true;
options.tstep_order = 1;
options.confined = true;

[options,prams] = initRigid2D(options,prams);

oc = curve;
xWalls = oc.createWalls(prams.Nbd, options);

% single fibre in centre of device
xc = [0;0];
tau = pi/2;

oc = curve;

Xfinal = rigid2D(options, prams, xc, tau, xWalls);

pp = post(['../output/data/',options.fileBase, '.dat'], ...
   ['../output/data/',options.fileBase,'_density.dat']);
pp.animated_gif('confined.gif', 1, [], 'fibres')
stats = pp.calculate_stats(1:prams.nv);
