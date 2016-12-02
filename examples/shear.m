close all

%% create animated gif and plots

prams.N = 128; % points per body
prams.Nbd = 0; %points on solid wall

g = -5;%shear rate

prams.nv = 8; % number of bodies
prams.nbd = 0; %number of walls
prams.lengths = ones(1,prams.nv);
prams.widths = 0.25*ones(1,prams.nv);
prams.T = 10*pi*(prams.lengths/prams.widths + prams.widths/prams.lengths)/abs(g);
prams.m = 800; % number of time steps
prams.order = 6;
prams.tracker_fnc = @(t) [20*cos(t),20*sin(t);5*cos(t),5*sin(t)];

options.farField = 'shear';
options.saveData = true;
options.fileBase = 'shear8';
options.append = false;
options.inear = true;
options.usePreco = false;
options.ifmm = true;
options.verbose = true;
options.profile = false;
options.tstep_order = 2;
options.confined = false;

[options,prams] = initRigid2D(options,prams);
xWalls = [];

xc = [linspace(0,10,prams.nv);zeros(1,prams.nv)];
tau  = pi/2*ones(1,prams.nv);

Xfinal = rigid2D(options, prams, xc, tau, xWalls);

pp = post(['../output/data/',options.fileBase, '.mat']);

a = prams.widths/2;
b = prams.lengths/2;
tau_exact = @(t) atan((b/a)*tan(a*b*g*t/(a^2+b^2)));
omega_exact = @(t) g/(a^2+b^2)*(b^2*cos(tau_exact(t)).^2 + a^2*sin(tau_exact(t)).^2);

pp.animated_gif('shear_8_fibers', 'tikz', 4, [], 'fibres')

