close all

prams.Np = 32; % points per body
prams.np = 2; % number of bodies
prams.T = 10; % time horizon
prams.number_steps = 200; % number of time steps
prams.lengths = 1;
prams.widths = 1;
prams.rounding_order = 2;
prams.minimum_separation = 1e-1;

options.far_field = 'extensional';
options.save = true;
options.file_base = 'extensional_test';
options.append = false;
options.near_singular = true;
options.use_precond = true;
options.fmm = true;
options.verbose = true;
options.profile = false;
options.tstep_order = 1;
options.confined = false;
options.resolve_collisions = true;
options.display_solution = true;
options.debug = false;
options.explicit = false;

[options,prams] = initRigid2D(options,prams);

tau = [pi/2, 0];
xc = [0,0; 2, 0];

Xfinal = rigid2DCollisions(options, prams, xc, tau, [], zeros(1,prams.np));