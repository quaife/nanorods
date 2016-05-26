prams.N = 32; % points per body
prams.nv = 49; % number of bodies
prams.T = 10; % time horizon
prams.m = 100; % number of time steps
prams.semimajors = 1*ones(1,prams.nv);
prams.semiminors = 0.25*ones(1,prams.nv);

options.farField = 'poiseuille';
options.usePlot = true;
options.axis = [-20 20 -5 5];
options.saveData = true;
options.dataFile = 'poiseuille_random_orientations';
options.append = false;
options.inear = true;
options.usePreco = true;

[options,prams] = initRigid2D(options,prams);

%% staggerd grid
x = linspace(0, 6*(sqrt(prams.nv)-1), 7);
y = linspace(0, 6*(sqrt(prams.nv)-1), 4);


[X1, Y1] = meshgrid(x,y);
[X2, Y2] = meshgrid(x + 6, y(2:end) - 6);

coeffr = 0;
xc = [[X1(:)', X2(:)'] + coeffr*(1 - 2*rand(1,prams.nv)); [Y1(:)'-7.5, Y2(:)'-7.5] + coeffr*(1 - 2*rand(1,prams.nv))];
tau = pi/2*ones(1,prams.nv) + pi*(1-2*rand(1,prams.nv));

%Xfinal = rigid2D(options, prams, xc, tau);

pp = post([options.dataFile,'.dat']);
pp.animated_gif('particles_staggered_49r_highv.gif', 1, [])
% stats = pp.calculate_stats(1:prams.nv);
% 
% theta = linspace(0,2*pi);
% pdf = zeros(length(stats.time),length(theta));
% 
% for i = 1:length(stats.time)    
%     for j = 1:length(theta)
%         pdf(i,j) = pp.recover_probability(theta(j), stats.a2(:,:,i), stats.a4(:,:,:,:,i));
%     end
% end
% 
% plot(theta(1:floor(end/2)), pdf(1,1:floor(end/2)));
% 
% hold
% plot(stats.hist_thetas(1:floor(end/2)), stats.pdf_actual(1,1:floor(end/2)), '-o');