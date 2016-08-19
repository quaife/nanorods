close all

prams.N = 32; % points per body
prams.Nbd = 192; %points on solid wall
prams.nv = 4; % number of bodies
prams.nbd = 2; %number of walls
prams.T = 20; % time horizon
prams.m = 200; % number of time steps
prams.lengths = 4.5*ones(1, prams.nv);
prams.widths = 1*ones(1,prams.nv);
prams.order = 2;
prams.tracker_fnc = @(t) [20*cos(t);20*sin(t);5*cos(4*t);-5*sin(4*t)];

options.farField = 'couette';
options.saveData = true;
options.fileBase = 'confined_flow';
options.append = false;
options.inear = true;
options.usePreco = false;
options.ifmm = true;
options.verbose = true;
options.profile = false;
options.tstep_order = 2;
options.confined = true;

[options,prams] = initRigid2D(options,prams);

oc = curve;
xWalls = oc.createWalls(prams.Nbd, options);

% single fibre in centre of device
xc = [9 0 -9 0; 0 11 0 -11];
tau = [0, pi/2, 0, pi/2];

oc = curve;

%Xfinal = rigid2D(options, prams, xc, tau, xWalls);

pp = post(['../output/data/',options.fileBase, '.dat'], ...
    ['../output/data/',options.fileBase, '_geometry','.dat'],...
    ['../output/data/',options.fileBase,'_density.dat']);
pp.animated_gif('confined_2nd_order_fibres.gif', 1, [], 'fibres')
%stats = pp.calculate_stats(1:prams.nv);
