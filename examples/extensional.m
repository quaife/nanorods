close all

prams.Np = 64; % points per body
prams.np = 2; % number of bodies
prams.T = 10; % time horizon
prams.number_steps = 200; % number of time steps
prams.lengths = 1;
prams.widths = 0.125;
prams.order = 2;

options.far_field = 'extensional';
options.save = true;
options.file_base = 'extensional';
options.append = false;
options.near_singular = true;
options.use_precond = false;
options.fmm = true;
options.verbose = true;
options.profile = false;
options.tstep_order = 2;
options.confined = false;
options.gmres_tol = 1e-6;

[options,prams] = initRigid2D(options,prams);

xc = [0, 0; 0.5627, 0];
tau = [pi/2,0];

Xfinal = rigid2D(options, prams, xc, tau);

pp = post(['../output/data/',options.file_base, '.mat']);

gif_options.file_name = 'extensional';
gif_options.file_type = 'gif';
gif_options.plot_fluid = false;
gif_options.xmin = 'auto:all';
gif_options.xmax = 'auto:all';
gif_options.ymin = 'auto:all';
gif_options.ymax = 'auto:all';
gif_options.axis = true;
gif_options.itmax = 'all';
gif_options.stride = 1;
gif_options.contour_field = [];
gif_options.velocity_quiver = false;

pp.animatedGif(gif_options);
