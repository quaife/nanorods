close all

prams.Np = 32; % points per body
prams.np = 48; % number of bodies
prams.T = 20; % time horizon
prams.number_steps = 40; % number of time steps
prams.lengths = 0.62;
prams.widths = 0.16;
prams.rounding_order = 2;
prams.minimum_separation = 5e-1;

options.far_field = 'taylor-green';
options.save = true;
options.file_base = 'taylor_green_test';
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

tau = zeros(1,prams.np);
nx = 4;
ny = 12;
x = linspace(-1,1,nx);
y = linspace(-1,1,ny);

[X,Y] = meshgrid(x,y);

xc = [X(:)'; Y(:)'];

Xfinal = rigid2DCollisions(options, prams, xc, tau, [], zeros(1,prams.np));