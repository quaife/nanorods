close all

prams.N = 32; % points per body
prams.Nbd = 0; %points on solid wall

g = -5;%shear rate

prams.nv = 8; % number of bodies
prams.nbd = 0; %number of walls
prams.lengths = ones(1,prams.nv);
prams.widths = 0.25*ones(1,prams.nv);
prams.T = 10*pi*(prams.lengths/prams.widths + prams.widths/prams.lengths)/abs(g);
prams.m = 200; % number of time steps
prams.order = 4;
prams.tracker_fnc = @(t) [20*cos(t),20*sin(t);5*cos(t),5*sin(t)];

options.farField = 'shear';
options.saveData = true;
options.fileBase = 'shear';
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

gif_options.file_name = 'shear';
gif_options.file_type = 'gif';
gif_options.plot_fluid = false;
gif_options.xmin = 'auto:frame';
gif_options.xmax = 'auto:frame';
gif_options.ymin = 'auto:all';
gif_options.ymax = 'auto:all';
gif_options.axis = true;
gif_options.itmax = prams.m;
gif_options.stride = 1;
    
pp.animated_gif(gif_options)

