close all

prams.N = 16; % points per body
prams.Nbd = 192; %points on solid wall

prams.nv = 2; % number of bodies
prams.nbd = 2; %number of walls
prams.T = 4*pi; % time horizon
prams.m = 1;%4*pi/0.1; % number of time steps
prams.lengths = 0.1*ones(1, prams.nv);
prams.widths = 0.1*ones(1,prams.nv);
prams.order = 2;
prams.tracker_fnc = @(t) [10,0;5*cos(t),5*sin(t)];
prams.gmresTol = 1e-8;

options.farField = 'couette';
options.saveData = true;
options.append = false;
options.inear = true;
options.usePreco = true;
options.ifmm = true;
options.verbose = true;
options.profile = false;
options.tstep_order = 2;
options.confined = true;
options.rk_tol = 1e-2;
options.dt_min = 1e-5;
options.rk_max_up = 1.5;
options.rk_max_down = 0.25;
options.rk_safety = 0.9;

[options,prams] = initRigid2D(options,prams);

oc = curve;
xWalls = oc.createWalls(prams.Nbd, options);

xc = [];
tau = [];

N = 6;
x = linspace(-10,10,N);
y = linspace(-10,10,N);
xc = zeros(2,0);
tau  = [];

buffer = 1.1*prams.lengths(1)/2;
[X,Y] = meshgrid(x,y);

dist = (x(2)-x(1))-prams.lengths(1);
max_disp = 0.9*dist;
X = X + max_disp*rand(size(X));
Y = Y + max_disp*rand(size(Y));

% remove particles outside domain
for m = 1:N
    for n = 1:N

       r = sqrt((X(m,n))^2+(Y(m,n))^2);
       if (r > 5 + buffer && r < 10 - buffer)
           xc(:,end+1) = [X(m,n);Y(m,n)];
           tau(end+1) = 0;
       end
    end 
end

%compute volume fraction
prams.nv = length(tau);
vol_frac = prams.nv*(prams.lengths(1)/2)^2/(10^2-5^2);
disp(['nv = ', num2str(prams.nv), ', volume fraction = ' num2str(vol_frac)]);

options.fileBase = ['couette_',num2str(prams.nv)];

Xfinal = rigid2D(options, prams, xc, tau, xWalls);
pp = post(['../output/data/',options.fileBase, '.mat']);

gif_options.file_name = 'couette_pressure_dense_small';
gif_options.file_type = 'tikz';
gif_options.plot_fluid = false;
gif_options.xmin = -10.5;
gif_options.xmax = 10.5;
gif_options.ymin = -10.5;
gif_options.ymax = 10.5;
gif_options.axis = false;
gif_options.itmax = 'all';
gif_options.stride = 1;
gif_options.contour_field = 'pressure';
gif_options.velocity_quiver = false;
pp.animatedGif(gif_options);

prams.lengths = 0.5*ones(1, prams.nv);
prams.widths = 0.5*ones(1,prams.nv);

N = 6;
xc = [];
tau = [];

x = linspace(-10,10,N);
y = linspace(-10,10,N);
xc = zeros(2,0);
tau  = [];

buffer = 1.1*prams.lengths(1)/2;
[X,Y] = meshgrid(x,y);

dist = (x(2)-x(1))-prams.lengths(1);
max_disp = 0.9*dist;
X = X + max_disp*rand(size(X));
Y = Y + max_disp*rand(size(Y));

% remove particles outside domain
for m = 1:N
    for n = 1:N

       r = sqrt((X(m,n))^2+(Y(m,n))^2);
       if (r > 5 + buffer && r < 10 - buffer)
           xc(:,end+1) = [X(m,n);Y(m,n)];
           tau(end+1) = 0;
       end
    end 
end

%compute volume fraction
prams.nv = length(tau);
vol_frac = prams.nv*(prams.lengths(1)/2)^2/(10^2-5^2);
disp(['nv = ', num2str(prams.nv), ', volume fraction = ' num2str(vol_frac)]);

options.fileBase = ['couette_',num2str(prams.nv)];

Xfinal = rigid2D(options, prams, xc, tau, xWalls);
pp = post(['../output/data/',options.fileBase, '.mat']);
gif_options.file_name = 'couette_pressure_dense_large';
pp.animatedGif(gif_options);