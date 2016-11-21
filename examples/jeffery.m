close all

%% create animated gif and plots

if 0
    prams.N = 64; % points per body
    prams.Nbd = 0; %points on solid wall

    g = -5;%shear rate

    prams.nv = 1; % number of bodies
    prams.nbd = 0; %number of walls
    prams.lengths = 1;
    prams.widths = 0.25;
    prams.T = pi*(prams.lengths/prams.widths + prams.widths/prams.lengths)/abs(g);
    prams.m = 400; % number of time steps
    prams.order = 2;
    prams.tracker_fnc = @(t) [20*cos(t),20*sin(t);5*cos(t),5*sin(t)];

    options.farField = 'shear';
    options.saveData = true;
    options.fileBase = 'jeffery';
    options.append = false;
    options.inear = true;
    options.usePreco = true;
    options.ifmm = true;
    options.verbose = true;
    options.profile = false;
    options.tstep_order = 2;
    options.confined = false;

    [options,prams] = initRigid2D(options,prams);
    xWalls = [];

    xc = [0;0];
    tau  = pi/2;

    %Xfinal = rigid2D(options, prams, xc, tau, xWalls);

    pp = post(['../output/data/',options.fileBase, '.mat']);

    a = prams.widths/2;
    b = prams.lengths/2;
    theta_exact = @(t) atan((b/a)*tan(a*b*g*t/(a^2+b^2)));
    omega_exact = @(t) g/(a^2+b^2)*(b^2*cos(theta_exact(t)).^2 + a^2*sin(theta_exact(t)).^2);

    % 
    % plot(pp.times, wrapTo2Pi(pp.tau));
    % hold on
    % plot(pp.times, pp.omega);
    % fplot(@(t) theta_exact(t) + pi/2, [0, prams.T]);
    % fplot(@(t) omega_exact(t), [0,prams.T]);
    % 
    % legend('angle', 'angular velocity', 'Jeffery angle', 'Jeffery velocity');


    %pp.animated_gif('jefery_orbit.gif', 3, [], 'fluid')
end

%% convergence study
if 1
    
    M = 5;
    errors = zeros(M,1);
    for i = 1:M
        
        disp(['dT = ', num2str(1/(2^(i+1)))]);
        prams.N = 256; % points per body
        prams.Nbd = 0; %points on solid wall

        g = -5;%shear rate

        prams.nv = 1; % number of bodies
        prams.nbd = 0; %number of walls
        prams.lengths = 1;
        prams.widths = 0.25;
        prams.T = 1;
        prams.m = 2^(i + 1); % number of time steps
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
        theta_exact = @(t) atan((b/a)*tan(a*b*g*t/(a^2+b^2)));
        omega_exact = @(t) g/(a^2+b^2)*(b^2*cos(theta_exact(t)).^2 + a^2*sin(theta_exact(t)).^2);
        
        errors(i) = abs(omega_exact(prams.T) - pp.omega(end));
    end
    
    rates = zeros(M-1,1);
    for i = 1:M-1
        rates(i) = log(errors(i)/errors(i+1))/log(2);
    end
end