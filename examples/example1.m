prams.N = 128; % points per body
prams.nv = 2; % number of bodies
prams.T = 10; % time horizon
prams.m = 100; % number of time steps

options.farField = 'extensional';
options.usePlot = true;
options.axis = [-10 10 -5 5];


[options,prams] = initRigid2D(options,prams);
% set options and parameters to default values

theta = (0:prams.N-1)'*2*pi/prams.N;
X = zeros(2*prams.N,prams.nv);
X(:,1) = [cos(theta) - 8;3*sin(theta)];
X(:,2) = [cos(theta) + 8;3*sin(theta)];
%X(:,3) = [cos(theta) - 8;3*sin(theta)];
%X(:,4) = [cos(theta) + 4;3*sin(theta)];
% set initial configuration

Xfinal = rigid2D(X,options,prams);

