close all
load ../output/data/couette_010_3_1.mat

prams.Np = 32; % points per body
prams.np = 160; % number of bodies
prams.Nw = prams.Np*8;
prams.nw = 2;

prams.T = 3*pi; % time horizon
prams.number_steps = 1000; % number of time steps
prams.lengths = 0.75;
prams.widths = 0.25;
prams.rounding_order = 2;
prams.minimum_separation = 1e-1;

options.far_field = 'couette';
options.save = true;
options.file_base = 'couette_010_3_1_tmp';
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
options.couette_speed = 2; % each revolution takes pi
options.explicit = false;
options.debug = false;

prams.tracker_fnc = @(t) [10*cos(2*t),10*sin(2*t);5,0];

[options,prams] = initRigid2D(options,prams);

oc = curve;
xWalls = oc.createWalls(prams.Nw, options);

restart = true;
animation = true;

if ~animation
    if restart
        xc = xc(:,:,end);
        tau = tau(end,:);
    else
        
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
            seed = 1000;
            rng(seed);
            
            total_np = prams.np;
            
            xc = [5;5];
            tau = 2*pi*rand(1,1);
            prams.np = 1;
            
            geom = capsules(prams, xc, tau);
            
            for i = 2:total_np
                
                max_attempts = 1000;
                attempt = 0;
                prams.np = i;
                
                while attempt < max_attempts
                    
                    add = true;
                    
                    r = r_in + (r_out-r_in)*rand(1,1);
                    theta = 2*pi*rand(1,1);
                    
                    xc_test = [xc, [r*cos(theta);r*sin(theta)]];
                    tau_test = [tau, 2*pi*rand(1,1)];
                    
                    geom_test = capsules(prams, xc_test, tau_test);
                    
                    
                    [near,~] = geom_test.getZone(geom_test,1);
                    
                    rx = sqrt(geom_test.X(1:end/2,end).^2+geom_test.X(end/2+1:end,end).^2);
                    
                    for j = 1:prams.np
                        
                        if (length(near.nearFibers{j}) > 0 || min(rx) < r_in + buffer || max(rx) > r_out - buffer)
                            add = false;
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
                    add = false;
                    disp('Failed to insert fibre, try increasing maximum number of iterations');
                end
                
                if add
                    xc = xc_test;
                    tau = tau_test;
                end
                
            end
        end
        
    end
    
    % xc = [6; -1.5];
    % tau = pi/4;
    
    %compute volume fraction
    prams.np = length(tau);
    vol_frac = prams.np*(prams.lengths/2*prams.widths/2)/(10^2-5^2);
    disp(['np = ', num2str(prams.np), ', volume fraction = ' num2str(vol_frac)]);    
    
    Xfinal = rigid2DCollisions(options, prams, xc, tau, xWalls);
end

pp = post(['../output/data/',options.file_base, '.mat']);

gif_options.file_name = 'couette_010_3_1';
gif_options.file_type = 'gif';
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
pp.animatedGif(gif_options);

