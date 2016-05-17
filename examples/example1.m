prams.N = 192; % points per body
prams.nv = 2; % number of bodies
prams.T = 20; % time horizon
prams.m = 100; % number of time steps
prams.semimajors = [3,3];
prams.semiminors = [1,1];

options.farField = 'shear';
options.usePlot = true;
options.axis = [-20 20 -5 5];
options.saveData = true;
options.dataFile = 'velocities.dat';
options.append = false;

xc = [-5, 5; 0.3, -0.3]; % [x-coordinates; y-coordinates]
tau = [pi/2, pi/2];

Xfinal = rigid2D(options, prams, xc, tau);

pp = post(options.dataFile);
pp.animated_gif('shear256.gif', 1)

