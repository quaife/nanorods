close all

prams.N = 64; % points per body
prams.nv = 2; % number of bodies
prams.T = 10; % time horizon
prams.m = 200; % number of time steps
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

xc = [0, 0.001; 0.6, -0.6];
tau = [pi/2,pi/2];

Xfinal = rigid2D(options, prams, xc, tau);

pp = post(['../output/data/',options.fileBase, '.mat']);

gif_options.file_name = 'extensional';
gif_options.file_type = 'gif';
gif_options.plot_fluid = false;
gif_options.xmin = 'auto:all';
gif_options.xmax = 'auto:all';
gif_options.ymin = 'auto:all';
gif_options.ymax = 'auto:all';
gif_options.axis = true;
gif_options.itmax = prams.m;
gif_options.stride = 1;
    
pp.animated_gif(gif_options);