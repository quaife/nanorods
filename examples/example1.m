%close all

prams.N = 48; % points per body
prams.nv = 25; % number of bodies
prams.T = 10; % time horizon
prams.m = 1000; % number of time steps
prams.lengths = 4.5*ones(1, prams.nv);
prams.widths = 1*ones(1,prams.nv);
prams.order = 4;

options.farField = 'extensional';
options.saveData = true;
options.fileBase = ['extensional_',datestr(now, 'yyyymmddHHMM')];
options.append = false;
options.inear = true;
options.usePreco = true;
options.ifmm = true;
options.verbose = true;

[options,prams] = initRigid2D(options,prams);

%% staggerd grid
% rown = 5;
% coln = prams.nv/(2*rown);
% 
% x = linspace(0, 3*max(prams.widths)*rown, rown);
% y = linspace(0, 1*max(prams.lengths)*coln, coln);
% 
% [X1, Y1] = meshgrid(x,y);
% [X2, Y2] = meshgrid(x - 2*max(prams.widths), y + max(prams.lengths)/2);
% 
% coeffr = 0;
% % xc = [[X1(:)', X2(:)'] + coeffr*(1 - 2*rand(1,prams.nv)); ...
% %     [Y1(:)' - coln/2*max(prams.lengths), Y2(:)'- coln/2*max(prams.lengths)] + coeffr*(1 - 2*rand(1,prams.nv))];
% xc = [[X1(:)', X2(:)'] + coeffr*(1 - 2*rand(1,prams.nv)); ...
%      [Y1(:)', Y2(:)'] + coeffr*(1 - 2*rand(1,prams.nv))];
% tau = pi/2*ones(1,prams.nv) + 2*coeffr*(1-2*rand(1,prams.nv));

%% non staggered grid, alternating orientations

rown = 5;
coln = prams.nv/rown;

x = linspace(-rown/2*max(prams.lengths), rown/2*max(prams.lengths), rown);
y = linspace(-coln/2*max(prams.lengths), coln/2*max(prams.lengths), coln);

coeffr = 1e-1;
[X, Y] = meshgrid(x,y);
xc = [X(:)' + coeffr*rand(1,prams.nv); Y(:)'+ (1 -2*coeffr*rand(1,prams.nv))];

tau = zeros(1,prams.nv);
tau(1:2:end) = pi/2;
    
Xfinal = rigid2D(options, prams, xc, tau);

pp = post([om.OUTPUTPATH_DATA,  options.dataFile, '.dat']);
pp.animated_gif('extenstional_large_timestep.gif', 1, [])
%stats = pp.calculate_stats(1:prams.nv);
