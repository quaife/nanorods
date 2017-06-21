close all

prams.Np = 32; % points per body
prams.np = 2; % number of bodies
prams.T = 5; % time horizon
prams.number_steps = 500; % number of time steps
prams.lengths = 1;
prams.widths = 0.125;
prams.rounding_order = 2;
prams.minimum_separation = 5e-1;

options.far_field = 'extensional';
options.save = true;
options.file_base = 'extensional_rk45';
options.append = false;
options.near_singular = true;
options.use_precond = true;
options.fmm = true;
options.verbose = true;
options.profile = false;
options.tstep_order = 2;
options.confined = false;
options.rk_tol = 1e-2;
options.rk_dt_min = 1e-2;
options.rk_max_up = 1.5;
options.rk_max_down = 0.25;
options.rk_safety = 0.9;
options.resolve_collisions = true;
options.display_solution = true;
options.debug = false;

[options,prams] = initRigid2D(options,prams);

% xc = [-1.95, -1.3, -0.65, 0, 0.65, 1.3, 1.95; 0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5];
% tau = [pi/2,pi/2,pi/2,pi/2,pi/2,pi/2,pi/2];
xc = [0,0 ; 0.6, 0];
tau = [0.95*pi/2, 0];
% xc = [-0.6, -0.2, 0.2, 0.6; 0, 0, 0, 0];
% tau = pi/4*ones(1,prams.np);
%xc = [1;0];
%tau = pi/2*ones(1,prams.np);

Xfinal = rigid2DCollisions(options, prams, xc, tau);

pp = post(['../output/data/',options.file_base, '.mat']);

gif_options.file_name = 'collision_resolved_N64_ellipse1';
gif_options.file_type = 'gif';
gif_options.plot_fluid = false;
gif_options.plot_pressure = false;
gif_options.xmin = 'auto:all';
gif_options.xmax = 'auto:all';
gif_options.ymin = 'auto:all';
gif_options.ymax = 'auto:all';
gif_options.axis = true;
gif_options.itmax = prams.number_steps;
gif_options.stride = 1;
gif_options.contour_field = [];    
gif_options.velocity_quiver = false;
pp.animatedGifExtended(gif_options);