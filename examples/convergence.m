rot = @(x, y, xc, yc) 1/((x-xc)^2+(y-yc)^2)*[(y-yc), -(x-xc)];

prams.N = 128; % points per body
prams.nv = 2; % number of bodies

prams.semimajors = [1,1];
prams.semiminors = [3,3];

options.farField = 'shear';
options.saveData = false;
options.dataFile = [];
options.append = false;
options.inear = true;

[options,prams] = initRigid2D(options,prams);

xc = [-1, 0; 0, 4]; % [x-coordinates; y-coordinates]
tau = [0, pi/2];
tt = tstep(options,prams);
om = monitor(options,prams);

%% single solve
geom = capsules(prams, xc, tau);
[density,Up,wp,iter,flag] = tt.timeStep(geom);


%% plot double layer potential

x = linspace(-10,10,30);
y = linspace(-10,10,30);

[X, Y] = meshgrid(x,y);
epsilon = 0.1;

om.plotDLP(geom, density, X,  Y, epsilon);
