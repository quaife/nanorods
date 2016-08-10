close all

prams.N = 32; % points per body
prams.Nbd = 192; %points on solid wall
prams.nv = 25; % number of bodies
prams.T = 1; % time horizon
prams.m = 10; % number of time steps
prams.lengths = 4.5*ones(1, prams.nv);
prams.widths = 1*ones(1,prams.nv);
prams.order = 4;

options.farField = 'circle';
options.saveData = true;
options.fileBase = 'confined_flow';
options.append = false;
options.inear = true;
options.usePreco = true;
options.ifmm = true;
options.verbose = true;
options.profile = true;
options.tstep_order = 1;
options.confined = true;

[options,prams] = initRigid2D(options,prams);

rng(123456); % set random seed


%% non staggered grid, alternating orientations

rown = 5;
coln = prams.nv/rown;

x = linspace(-rown/2*max(prams.lengths), rown/2*max(prams.lengths), rown);
y = linspace(-coln/2*max(prams.lengths), coln/2*max(prams.lengths), coln);

coeffr = 1e-1;
[X, Y] = meshgrid(x,y);
xc = [X(:)' + coeffr*rand(1,prams.nv); Y(:)'+ (1 -2*coeffr*rand(1,prams.nv))];

tau = zeros(1,prams.nv);
tau(1:2:end) = pi/2;

oc = curve;
xWalls = oc.createWalls(prams.Nbnd, options);

Xfinal = rigid2D(options, prams, xc, tau, xWalls);

 pp = post(['../output/data/',options.fileBase, '.dat'], ...
   ['../output/data/',options.fileBase,'_density.dat']);
pp.animated_gif('confined.gif', 1, [], 'fibres')
%pp.animated_gif('extenstional_fluid_second_order.gif', 1, [], 'fluid')
stats = pp.calculate_stats(1:prams.nv);
