function [] = couette(base_name, np, lengths, widths, rotation_speed, n_revolutions, restart)

%load(['../output/data/', fname, '.mat']);

%MODIFY THESE PARAMETERS
prams.np = np; % number of bodies
prams.lengths = lengths;
prams.widths = widths;
options.couette_speed = rotation_speed;
prams.T = n_revolutions*2*pi/options.couette_speed;  
%prams.tracker_fnc = [];
prams.tracker_fnc = @(t) [10*cos(1*t),10*sin(1*t);5,0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prams.Np = 32; % points per body
prams.Nw = prams.Np*8;
prams.nw = 2;

prams.number_steps = n_revolutions*500; % number of time steps
prams.rounding_order = 2;
prams.minimum_separation = 1e-1;

options.far_field = 'couette';
options.save = true;
options.append = false;
options.near_singular = true;
options.use_precond = true;
options.fmm = true;
options.verbose = true;
options.profile = false;
options.tstep_order = 1;
options.confined = true;
options.resolve_collisions = true;
options.display_solution = true;
options.explicit = false;
options.debug = false;
options.sedementation = false;

[options,prams] = initRigid2D(options,prams);
options.file_base = 'couette_einstein_80';
disp(['COUETTE SPEED: ', num2str(options.couette_speed)]);

options.file_base = base_name;
 

oc = curve;
animation = true;

if ~animation
    if restart
        load(['../output/data/',base_name,'.mat']);
        prams.minimum_separation = 1e-1;
        
        prams.Np = 32; % points per body
        prams.Nw = prams.Np*8;

        xWalls = oc.createWalls(prams.Nw, options);
        xc = xc(:,:,1:end);
        tau = tau(1:end,:);

    else
        xWalls = oc.createWalls(prams.Nw, options);
        xc = zeros(2,0);
        tau  = [];
        r_in = 5;
        r_out = 10;
        buffer = 0.5*prams.lengths;
        
        structured = false;
        
        if structured % ADD FIBRES ON REGULAR GRID
            
            N = 3; % number of rings
            N_inner = 50;
            
            buffer = 0.75*prams.lengths;
            r = linspace(r_in + buffer, r_out - buffer, N);
            
            spacing = ((r_out - r_in) - prams.lengths*N)/(N+1);
            
            for i = 1:length(r)
                
                N = round(2*pi*(r(i)-r_in)+N_inner);
                R = r(i)*ones(1,N);
                
                start = mod(i,2)  + 1;
                
                n = randi([-1,1], 1, N);
                
                R(start:2:end) = R(start:2:end) + 0.4*n(start:2:end)*spacing;
                
                THETA = linspace(0,2*pi - 2*pi/N,N);
                
                %[R,THETA] = meshgrid(radius, theta);
                
                % R = R + 0.8*spacing*rand(size(R));
                X = R.*cos(THETA);
                Y = R.*sin(THETA);
                
                for m = 1:size(X,1)
                    for n = 1:size(X,2)
                        xc(:,end+1) = [X(m,n);Y(m,n)];
                        tau(end+1) = atan2(Y(m,n),X(m,n));
                    end
                end
            end
            
        else % ADD FIBERS BY MONTE CARLO
            seed = 1001;
            rng(seed);
            
            total_np = prams.np;
            walls = capsules([], xWalls);
            prams_test = prams;
            prams_test.np = 1;
            
            for i = 1:total_np
                
                max_attempts = 5000;
                attempt = 0;
                prams.np = i;
                
                while attempt < max_attempts
                    
                    add = true;
                    
                    r = r_in + (r_out-r_in)*rand(1,1);
                    theta = 2*pi*rand(1,1);
                    
                    xc_test = [r*cos(theta);r*sin(theta)];
                    tau_test = 2*pi*rand(1,1);
                    
                    particle_test = capsules(prams_test, xc_test, tau_test);
                    [~, nearw] = walls.getZone(particle_test,2);
                    
                    rx = sqrt(particle_test.X(1:end/2,end).^2+particle_test.X(end/2+1:end,end).^2);
                    
                    if max(rx) > r_out || min(rx) < r_in
                        add = false;
                    else
                        
                        if i > 1
                            
                            [~, nearp] = particle_test.getZone(geom,2);
                            
                            if length(nearp.nearFibers{1}) > 0 || length(vertcat(nearw.nearFibers{:})) > 0
                                add = false;
                            end
                        else
                            
                            if  length(nearw.nearFibers{1}) > 0
                                add = false;
                            end
                        end
                        
                    end
                    
                    attempt = attempt + 1;
                    
                    if add
                        disp(['inserted fibre ', num2str(i)]);
                        break;
                    else
                        disp(['failed to insert fibre ', num2str(i), '(attempt ', num2str(attempt), '), trying again']);
                    end
                end
                
                if (attempt == max_attempts)
                    disp('Failed to insert desired number of fibres, try increasing maximum number of iterations');
                    break;
                end
                
                if add
                    xc = [xc, xc_test];
                    tau = [tau, tau_test];
                    
                   geom = capsules(prams, xc, tau);
                end
                
            end
        end
        
    end
    
    % xc = [6; -1.5];
    % tau = pi/4;
    
    %compute volume fraction
    prams.np = size(tau,2);
    vol_frac = prams.np*(prams.lengths/2*prams.widths/2)/(10^2-5^2);
    disp(['np = ', num2str(prams.np), ', volume fraction = ' num2str(vol_frac)]);    
    
    Xfinal = rigid2DCollisions(options, prams, xc, tau, xWalls, zeros(1,prams.np));
end

if animation
    pp = post(['../output/data/',options.file_base, '.mat']);

    pp.prams.tracker_fnc = @(t) [10*cos(1*t),10*sin(1*t);5,0];
    gif_options.file_name = options.file_base;
    gif_options.file_type = 'tikz';
    gif_options.plot_fluid = false;
    gif_options.xmin = -10.5;
    gif_options.xmax = 10.5;
    gif_options.ymin = -10.5;
    gif_options.ymax = 10.5;
    gif_options.axis = false;
    gif_options.itmax = 'all';
    gif_options.stride = 1;
    gif_options.contour_field = [];
    gif_options.velocity_quiver = false;
    gif_options.dt = [];
    gif_options.tracers = false;
    pp.animatedGif(gif_options);
 end
