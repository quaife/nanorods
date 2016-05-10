prams.N = 32; % points per body
prams.nv = 2; % number of bodies
prams.T = 1; % time horizon
prams.m = 100; % number of time steps
prams.gmresTol = 1e-8; % GMRES tolerance

options.order = 1; 
% time stepping order (only options is currently first)
options.inear = false; % flag for near-singular integration

addpath ../src
% TODO: Need a routine to set the undefined prams and options

theta = (0:prams.N-1)'*2*pi/prams.N;
X = zeros(2*prams.N,prams.nv);
X(:,1) = [cos(theta);3*sin(theta)];
X(:,2) = [cos(theta) - 4;3*sin(theta)+3];
% set initial configuration

tt = tstep(options,prams);

[Xnew,iter,iflag] = tt.timeStep(X);

