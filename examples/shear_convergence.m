close all

% create fine solution
prams.Nw = 0; %points on solid wall

%g = -5;%shear rate

prams.np = 8; % number of bodies
prams.nw = 0; %number of walls
prams.lengths = 1;
prams.widths = 0.25;
prams.T = 15;
prams.number_steps = 300; % number of time steps
prams.rounding_order = 2;
prams.tracker_fnc = @(t) [20*cos(t),20*sin(t);5*cos(t),5*sin(t)];
prams.minimum_separation = 0.1;

options.far_field = 'shear';
options.save = true;
options.file_base = 'shear_fine';
options.append = false;
options.near_singular = true;
options.use_precond = true;
options.fmm = true;
options.verbose = true;
options.profile = false;
options.tstep_order = 2;
options.confined = false;
options.resolve_collisions = false;
options.debug = false;
options.display_solution = true;

[options,prams] = initRigid2D(options,prams);
xWalls = [];

xc = [linspace(0,10,prams.np);zeros(1,prams.np)];
%xc = [0;0];
tau  = pi/2*ones(1,prams.np);

% create fine solution
prams.Np = 512; % points per body
tic;
Xfinal = rigid2DCollisions(options, prams, xc, tau, xWalls);
t_fine = toc;

pp = post(['../output/data/',options.file_base, '.mat']);
xc_fine = pp.xc;
tau_fine = pp.tau;

N = [8,16,32,64,128,256];
xc_err = zeros(1,length(N));
tau_err = zeros(1,length(N));
times = zeros(1,length(N));

for i = 1:length(N)
    prams.Np = N(i);
    options.file_base = ['shear_',num2str(N(i))];
    tic;
    Xfinal = rigid2DCollisions(options, prams, xc, tau, xWalls);
    times(i) = toc;
    
    pp = post(['../output/data/',options.file_base, '.mat']);
    xc_err(i) = norm(pp.xc(:,:,end) - xc_fine(:,:,end), 'fro');
    tau_err(i) = norm(pp.tau(end,:) - tau_fine(end,:));
end