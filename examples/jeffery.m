close all

%% create animated gif and plots
if 1
    prams.N = 64; % points per body
    prams.Nbd = 0; %points on solid wall

    g = -5;%shear rate

    prams.nv = 1; % number of bodies
    prams.nbd = 0; %number of walls
    prams.lengths = 1;
    prams.widths = 0.25;
    prams.T = pi*(prams.lengths/prams.widths + prams.widths/prams.lengths)/abs(g);
    prams.m = 100; % number of time steps
    prams.order = 2;
    prams.tracker_fnc = @(t) [20*cos(t),20*sin(t);5*cos(t),5*sin(t)];

    options.farField = 'shear';
    options.saveData = true;
    options.fileBase = 'jeffery';
    options.append = false;
    options.inear = true;
    options.usePreco = false;
    options.ifmm = true;
    options.verbose = true;
    options.profile = false;
    options.tstep_order = 2;
    options.confined = false;

    [options,prams] = initRigid2D(options,prams);
    xWalls = [];

    xc = [0;0];
    tau  = pi/2;

    Xfinal = rigid2D(options, prams, xc, tau, xWalls);

    pp = post(['../output/data/',options.fileBase, '.mat']);

    a = prams.widths/2;
    b = prams.lengths/2;
    tau_exact = @(t) atan((b/a)*tan(a*b*g*t/(a^2+b^2)));
    omega_exact = @(t) g/(a^2+b^2)*(b^2*cos(tau_exact(t)).^2 + a^2*sin(tau_exact(t)).^2);

    
    fplot(@(t) tau_exact(t) + pi/2, [0, prams.T]);
    hold on
    fplot(@(t) omega_exact(t), [0,prams.T]);
    plot(pp.times, wrapTo2Pi(pp.tau), 'bo');
    plot(pp.times(1:end-1), pp.omega, 'ro');
    
    legend({'$\theta$ (Jeffery)', '$\omega$ (Jeffery)', '$\theta$ (BIE)', '$\omega$ (BIE)'}, 'interpreter', 'latex');
    xlabel('$t$', 'interpreter', 'latex');

    addpath('../tests/matlab2tikz/src');
    matlab2tikz('jeffery.tex', 'height', '10cm', 'width', '12cm');
    
    pp.animated_gif('jeffery_orbit', 'tikz', 1, [], 'fluid')
end

%% convergence study
if 0
    
    M = 4;
    errors_fe = zeros(M,1);
    errors_ab = zeros(M,1);
    dt = zeros(M,1);
    
    for i = 1:M
        
        dt(i) = 1/(2^(i+3));
        disp(['dT = ', num2str(1/(2^(i+2)))]);
        prams.N = 128; % points per body
        prams.Nbd = 0; %points on solid wall

        g = -5;%shear rate

        prams.nv = 1; % number of bodies
        prams.nbd = 0; %number of walls
        prams.lengths = 1;
        prams.widths = 0.25;
        prams.T = 0.125;
        prams.m = prams.T*(2^(i + 3)); % number of time steps
        prams.order = 2;
        prams.tracker_fnc = @(t) [20*cos(t),20*sin(t);5*cos(t),5*sin(t)];

        options.farField = 'shear';
        options.saveData = true;
        options.fileBase = 'jeffery_convergence';
        options.append = false;
        options.inear = true;
        options.usePreco = false;
        options.ifmm = true;
        options.verbose = true;
        options.profile = false;
        options.tstep_order = 1;
        options.confined = false;

        [options,prams] = initRigid2D(options,prams);
        xWalls = [];

        xc = [0;0];
        tau  = pi/2;

        Xfinal = rigid2D(options, prams, xc, tau, xWalls);

        pp = post(['../output/data/',options.fileBase, '.mat']);

        a = prams.widths/2;
        b = prams.lengths/2;
        tau_exact = @(t) atan((b/a)*tan(a*b*g*t/(a^2+b^2)));
        omega_exact = @(t) g/(a^2+b^2)*(b^2*cos(tau_exact(t)).^2 + a^2*sin(tau_exact(t)).^2);
        
        errors_fe(i) = abs(tau_exact(pp.times(end)) - pp.tau(end) + pi/2);

        options.tstep_order = 2;
        Xfinal = rigid2D(options, prams, xc, tau, xWalls);

        pp = post(['../output/data/',options.fileBase, '.mat']);
        errors_ab(i) = abs(tau_exact(pp.times(end)) - pp.tau(end) + pi/2);
    end
    
    figure();
    loglog(dt, errors_fe, 'bo', 'linewidth', 2);
    hold on
    loglog(dt, errors_ab, 'ro', 'linewidth', 2);

    loglog([dt(1), dt(end)], [errors_fe(1), errors_fe(1)/(2^(length(dt)-1))], 'b');
    loglog([dt(1), dt(end)], [errors_ab(1), errors_ab(1)/(4^(length(dt)-1))], 'r');
    
    legend({'forward Euler', 'Adams Bashforth', '$O(\Delta t)$', '$O(\Delta t^2)$'}, 'interpreter', 'latex');
    xlabel('$\Delta t$', 'interpreter', 'latex');
    ylabel('error in $\theta$', 'interpreter', 'latex');
end
