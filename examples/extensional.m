close all

prams.N = 64; % points per body
prams.nv = 2; % number of bodies
prams.T = 5; % time horizon
prams.m = 100; % number of time steps
prams.lengths = 1*ones(1, prams.nv);
prams.widths = 1*ones(1,prams.nv);
prams.order = 2;

options.farField = 'extensional';
options.saveData = true;
options.fileBase = 'extensional';
options.append = false;
options.inear = true;
options.usePreco = false;
options.ifmm = true;
options.verbose = true;
options.profile = false;
options.tstep_order = 2;
options.confined = false;

[options,prams] = initRigid2D(options,prams);

xc = [0, 0.05; 0.6, -0.6];
tau = [pi/2,pi/2];

Xfinal = rigid2D(options, prams, xc, tau);

pp = post(['../output/data/',options.fileBase, '.mat']);
pp.animated_gif('extensional', 'gif', 1, [], 'fibres');