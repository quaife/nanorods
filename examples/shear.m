%close all
profile off

prams.N = 96; % points per body
prams.nv = 50; % number of bodies
prams.T = 10; % time horizon
prams.m = 100; % number of time steps
prams.lengths = 4.5*ones(1, prams.nv);
prams.widths = 1*ones(1,prams.nv);
prams.order = 4;

options.farField = 'shear';
options.saveData = true;
options.fileBase = 'shear_periodic';
options.append = false;
options.inear = true;
options.usePreco = true;
options.ifmm = true;
options.verbose = true;
options.profile = true;
options.tstep_order = 2;
options.n_cores_matlab = 2;

[options,prams] = initRigid2D(options,prams);

rng(123456); % set random seed

%% staggerd grid
rown = 5;
coln = 5;

x1 = linspace(-1.1*max(prams.lengths)*(rown - 1)/2, 1.1*max(prams.lengths)*(rown - 1)/2, rown);
y1 = x1;

[X1, Y1] = meshgrid(x1,y1);
X2 = X1 + (x1(2) - x1(1))/2;
Y2 = Y1 - (y1(2) - y1(1))/2;

coeffr = 0;
xc = [[X1(:)', X2(:)'] + coeffr*(1 - 2*rand(1,prams.nv)); ...
    [Y1(:)', Y2(:)'] + coeffr*(1 - 2*rand(1,prams.nv))];
tau = pi/2*ones(1,prams.nv) + 2*coeffr*(1-2*rand(1,prams.nv));
    
Xfinal = rigid2DPeriodic(options, prams, xc, tau, min(xc(1,:)), max(xc(1,:)), min(xc(2,:)), max(xc(2,:)), (x1(2)-x1(1))/2);

pp = post(['../output/data/',options.fileBase, '.dat'], ...
    ['../output/data/',options.fileBase,'_density.dat']);
% 
% for k = 7:-1:1
%pp.plot_fibres(105, -15, 15, -15, 15)
% pause
% end
pp.animated_gif('testVerySmall.gif', 10, [], 'fibres')
%pp.animated_gif('extenstional_fluid_second_order.gif', 1, [], 'fluid')
%stats = pp.calculate_stats(1:prams.nv);

