close all

%% create animated gif and plots
if 0
    prams.Np = 64; % points per body
    prams.Nw = 0; %points on solid wall

    g = 1;%shear rate

    prams.np = 1; % number of bodies
    prams.nw = 0; %number of walls
    prams.lengths = 3;
    prams.widths = 1;
    prams.T = 1*pi/2*(prams.lengths/prams.widths + prams.widths/prams.lengths)/abs(g);
    prams.number_steps = 50; % number of time steps
    prams.order = 2;

    options.far_field = 'shear';
    options.save_data = true;
    options.file_base = 'jeffery';
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
    options.explicit = false;

    [options,prams] = initRigid2D(options,prams);
    xWalls = [];

    xc = [0;0];
    tau  = pi/2;

    Xfinal = rigid2DCollisions(options, prams, xc, tau, xWalls, 0);

    pp = post(['../output/data/',options.file_base, '.mat']);

    if 1
        % make animated gif
        gif_options.file_name = 'jeffery_orbit_dissipation';
        gif_options.file_type = 'tikz';
        gif_options.velocityQuiver = false;
        gif_options.xmin = -2;
        gif_options.xmax = 2;
        gif_options.ymin = -1;
        gif_options.ymax = 1;
        gif_options.axis = false;
        gif_options.itmax = prams.m;
        gif_options.stride = 1;
        gif_options.grid_pts = 50;
        gif_options.velocity_grid_stride = 1;
        gif_options.bg_flow = @(x,y) [y, zeros(length(y),1)];
        gif_options.contour_field = 'dissipation';
        
        % pp.animatedGif(gif_options);
        
        % compare with exact solution
        a = prams.widths/2;
        b = prams.lengths/2;
        tau_exact = @(t) atan((b/a)*tan(a*b*g*t/(a^2+b^2)));
        omega_exact = @(t) g/(a^2+b^2)*(b^2*cos(tau_exact(t)).^2 + a^2*sin(tau_exact(t)).^2);
        
        figure();
        fplot(@(t) tau_exact(t) + pi/2, [0, prams.T]);
        hold on
        fplot(@(t) omega_exact(t), [0,prams.T]);
        plot(pp.times, wrapTo2Pi(pp.tau), 'bo');
        plot(pp.times(1:end-1), pp.omega, 'ro');
        
        legend({'$\theta$ (Jeffery)', '$\omega$ (Jeffery)', '$\theta$ (BIE)', '$\omega$ (BIE)'}, 'interpreter', 'latex');
        xlabel('$t$', 'interpreter', 'latex');
        
        addpath('../src/matlab2tikz/src');
        matlab2tikz('jeffery.tex', 'height', '10cm', 'width', '12cm');
    end
end

%% convergence study
if 1
    
    M = 4;
    errors_fe = zeros(M,1);
    errors_ab = zeros(M,1);
    dt = zeros(M,1);
    
    for i = 1:M
        
        dt(i) = 1/(2^(i+3));
        disp(['dT = ', num2str(1/(2^(i+2)))]);
        prams.Np = 128; % points per body
        prams.Nw = 0; %points on solid wall

        g = 1;%shear rate

        prams.np = 1; % number of bodies
        prams.nw = 0; %number of walls
        prams.lengths = 1;
        prams.widths = 0.25;
        prams.T = pi/4*(prams.lengths/prams.widths + prams.widths/prams.lengths)/abs(g);
        prams.number_steps = prams.T/dt(i); % number of time steps
        prams.order = 2;
        prams.minimum_separation = 0.1;
        
        options.far_field = 'shear';
        options.save_data = true;
        options.file_base = 'jeffery';
        options.append = false;
        options.near_singular = true;
        options.use_precond = true;
        options.fmm = true;
        options.verbose = true;
        options.profile = false;
        options.tstep_order = 1;
        options.confined = false;
        options.resolve_collisions = false;
        options.debug = false;
        options.display_solution = true;
        options.explicit = true;

        [options,prams] = initRigid2D(options,prams);
        xWalls = [];

        xc = [0;0];
        tau  = pi/2;

        Xfinal = rigid2DCollisions(options, prams, xc, tau, xWalls, 0);

        pp = post(['../output/data/',options.file_base, '.mat']);

        a = prams.widths/2;
        b = prams.lengths/2;
        tau_exact = @(t) atan((b/a)*tan(a*b*g*t/(a^2+b^2)));
        omega_exact = @(t) g/(a^2+b^2)*(b^2*cos(tau_exact(t)).^2 + a^2*sin(tau_exact(t)).^2);
        
        errors_fe(i) = abs(tau_exact(pp.times(end)) - pp.tau(end) + pi/2);

        options.tstep_order = 2;
        Xfinal = rigid2DCollisions(options, prams, xc, tau, xWalls,0);

        pp = post(['../output/data/',options.file_base, '.mat']);
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
