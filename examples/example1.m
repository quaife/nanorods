prams.N = 32; % points per body
prams.nv = 49; % number of bodies
prams.T = 10; % time horizon
prams.m = 100; % number of time steps
prams.semimajors = 2*ones(1,prams.nv);
prams.semiminors = 0.5*ones(1,prams.nv);

options.farField = 'shear';
options.usePlot = true;
options.axis = [-20 20 -5 5];
options.saveData = true;
options.dataFile = 'velocities.dat';
options.append = false;
options.inear = true;

[options,prams] = initRigid2D(options,prams);

xc = [linspace(0,5*prams.nv,prams.nv); ...
       rand(1,1)*1.1*ones(1,prams.nv/2), -rand(1,1)*1.1*ones(1,prams.nv/2)]; % [x-coordinates; y-coordinates]
%xc = [0;-1.1];
tau = 0*ones(1,prams.nv);

%Xfinal = rigid2D(options, prams, xc, tau);

pp = post(options.dataFile);
pp.animated_gif('shear10_particles64.gif', 1, [])
stats = pp.calculate_stats(1:prams.nv);
