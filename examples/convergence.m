rot = @(x, y, xc, yc) 1/((x-xc)^2+(y-yc)^2)*[(y-yc), -(x-xc)];

close all

prams.N = 16; % points per body
prams.nv = 2; % number of bodies

prams.semimajors = [1, 1];
prams.semiminors = [1, 1];

options.farField = 'rotlet';
options.saveData = false;
options.dataFile = [];
options.append = false;
options.inear = true;

[options,prams] = initRigid2D(options,prams);

xc = [0 4;0 4]; % [x-coordinates; y-coordinates]
tau = [pi/4 0];
tt = tstep(options,prams);
om = monitor(options,prams);

%% single solve
geom = capsules(prams, xc, tau);
[density,Up,wp,iter,flag] = tt.timeStep(geom);


%% plot double layer potential

x = linspace(-10,10,30);
y = linspace(-10,10,30);

[X, Y] = meshgrid(x,y);
epsilon = 0.5;

U = om.plotDLP(geom, density, X,  Y, epsilon);
