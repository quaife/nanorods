close all

prams.N = 32; % points per body
prams.nv = 49; % number of bodies
prams.T = 10; % time horizon
prams.m = 100; % number of time steps
prams.semimajors = 2*ones(1,prams.nv);
prams.semiminors = 0.5*ones(1,prams.nv);

options.farField = 'poiseuille';
options.usePlot = true;
options.axis = [-20 20 -5 5];
options.saveData = true;
<<<<<<< HEAD
options.dataFile = 'velocities_49_particles_staggeredr_highv';
=======
%options.dataFile = 'velocities_49_particles_staggeredr';
>>>>>>> 34af72a5f0ee67ad239b66978662cb71636053fe
options.append = false;
options.inear = true;
options.usePreco = false;

[options,prams] = initRigid2D(options,prams);


% x = linspace(-1*sqrt(prams.nv), 1*sqrt(prams.nv), sqrt(prams.nv));
% y = linspace(-2*sqrt(prams.nv), 2*sqrt(prams.nv), sqrt(prams.nv));
% 
% [X, Y] = meshgrid(x,y);
% 
% xc = [X(:)'+ 1*rand(1,prams.nv); Y(:)'+0.5*rand(1,prams.nv)];
% tau = pi/2*ones(1,prams.nv);

%% staggerd grid
x = linspace(0, 3*(sqrt(prams.nv)-1), 7);
y = linspace(0, 2.5*(sqrt(prams.nv)-1), 4);


[X1, Y1] = meshgrid(x,y);
[X2, Y2] = meshgrid(x + 1.5, y(2:end) - 2.5);

coeffr = 0.1;
xc = [[X1(:)', X2(:)'] + coeffr*(1 - 2*rand(1,prams.nv)); [Y1(:)'-7.5, Y2(:)'-7.5] + coeffr*(1 - 2*rand(1,prams.nv))];
tau = pi/2*ones(1,prams.nv) + 2*coeffr*(1-2*rand(1,prams.nv));

% xc = [linspace(0,5*prams.nv,prams.nv); ...
%        rand(1,1)*1.1*ones(1,prams.nv/2), -rand(1,1)*1.1*ones(1,prams.nv/2)]; % [x-coordinates; y-coordinates]
%xc = [0;-1.1];
%tau = 0*ones(1,prams.nv);

Xfinal = rigid2D(options, prams, xc, tau);

pp = post([options.dataFile,'.dat']);
pp.animated_gif('particles_staggered_49r_highv.gif', 1, [])
stats = pp.calculate_stats(1:prams.nv);
