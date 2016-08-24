close all

prams.N = 48; % points per body
prams.Nbd = 192; %points on solid wall
prams.nv = 4; % number of bodies
prams.nbd = 2; %number of walls
prams.T = 20; % time horizon
prams.m = 400; % number of time steps
prams.lengths = 1*ones(1, prams.nv);
prams.widths = 0.5*ones(1,prams.nv);
prams.order = 2;
prams.tracker_fnc = @(t) [20*cos(t),20*sin(t);5*cos(t),5*sin(t)];
%prams.tracker_fnc = @(t) [10*cos(t), 10*sin(t);5*cos(t), 5*sin(t)];

options.farField = 'couette';
options.saveData = true;
options.fileBase = 'couette';
options.append = false;
options.inear = true;
options.usePreco = true;
options.ifmm = true;
options.verbose = true;
options.profile = false;
options.tstep_order = 2;
options.confined = true;

[options,prams] = initRigid2D(options,prams);

oc = curve;
xWalls = oc.createWalls(prams.Nbd, options);

xc = [9 0 -9 0; 0 7 0 -7];
tau = [0, pi/2, 0, pi/2];
% tau = 2*pi*rand(1,1);
% xc = [6;6];

%Xfinal = rigid2D(options, prams, xc, tau, xWalls);

pp = post(['../output/data/',options.fileBase, '.mat']);
pp.animated_gif('couette.gif', 1, [], 'fibres')
%stats = pp.calculate_stats(1:prams.nv);
