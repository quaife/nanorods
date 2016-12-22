close all

prams.N = 32; % points per body
prams.Nbd = 192; %points on solid wall

prams.nv = 2; % number of bodies
prams.nbd = 2; %number of walls
prams.T = 10; % time horizon
prams.m = 100; % number of time steps
prams.lengths = 0.5*ones(1, prams.nv);
prams.widths = 0.5*ones(1,prams.nv);
prams.order = 2;
prams.tracker_fnc = @(t) [20*cos(t),20*sin(t);5*cos(t),5*sin(t)];
prams.gmresTol = 1e-8;
prams.a = [0, 0, 0, 0, 0;...
           1/4, 0, 0, 0, 0;...
           3/32, 9/32, 0, 0 ,0;...
           1932/2197, -7200/2197, 7296/2197, 0, 0;...
           439/216, -8, 3680/513, -845/4104, 0;
           -8/27, 2, -3544/2565, 1859/4104, -11/40];
prams.b1 = [25/216;0;1408/2565;2197/4104;-1/5];
prams.b2 = [16/135;0;6656/12825;28561/56430;-9/50;2/55];

prams.tol = 1e-5;
prams.dt_min = 1e-5;

options.farField = 'couette';
options.saveData = true;
options.fileBase = 'couette';
options.append = false;
options.inear = true;
options.usePreco = true;
options.ifmm = true;
options.verbose = true;
options.profile = false;
options.tstep_order = 'adaptive';
options.confined = true;

[options,prams] = initRigid2D(options,prams);

oc = curve;
xWalls = oc.createWalls(prams.Nbd, options);

% add particles on uniform grid
N = 32;
x = linspace(-10,10,N);
y = linspace(-10,10,N);

buffer = 0.3;
[X,Y] = meshgrid(x,y);
xc = zeros(2,0);
tau  = [];

% remove particles outside domain
for i = 1:N
    for j = 1:N
        
       r = sqrt(X(i,j)^2+Y(i,j)^2);
       if (r > 5 + buffer && r < 10 - buffer)
           xc(:,end+1) = [X(i,j);Y(i,j)];
           tau(end+1) = 0;
       end
    end 
end

% compute volume fraction
prams.nv = length(tau);
vol_frac = prams.nv*(prams.lengths(1)/2)^2/(10^2-5^2);
disp(['nv = ', num2str(prams.nv), ', volume fraction = ' num2str(vol_frac)]);


Xfinal = rigid2DAdaptive(options, prams, xc, tau, xWalls);

pp = post(['../output/data/',options.fileBase, '.mat']);

gif_options.file_name = 'couette';
gif_options.file_type = 'gif';
gif_options.plot_fluid = false;
gif_options.xmin = -10.5;
gif_options.xmax = 10.5;
gif_options.ymin = -10.5;
gif_options.ymax = 10.5;
gif_options.axis = false;
gif_options.itmax = prams.m;
gif_options.stride = 1;
gif_options.plot_pressure = false;
gif_options.grid_pts = 400;

pp.animated_gif(gif_options);

