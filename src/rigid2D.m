function [Xfinal] = rigid2D(options, prams, xc, tau, xWalls)

% delete(gcp('nocreate'))
% parpool(options.n_cores_matlab)   

ttotal = tic;
    

geom = capsules(prams, xc, tau);

if (options.confined)
    walls = capsules(prams, xWalls);
else
    walls = [];
end

% [xc_added, tau_added] = geom.fill_couette(5, 10, 86, prams);
% xc = [xc, xc_added];
% tau = [tau, tau_added];

% prams.nv = length(tau);
% prams.legnths = prams.lengths(1)*ones(1,prams.nv);
% prams.widths = prams.widths(1)*ones(1,prams.nv);

om = monitor(options, prams, xc, tau);
tt = tstep(options, prams, om, geom, walls, tau);
potF = poten(geom.N,om);

if (om.profile)
    profile on;
end

% read in existing file if necessary
if options.append
    
    pp = post('../output/data/',options.fileBase, '.mat');

    time = pp.times(end);
    iT = length(pp.times);
    
    xc = pp.xc(end,:);
    tau = pp.tau(end,:);    
    
    if (options.tstep_order == 2)
        Up_m1 = pp.U(end,:);
        wp_m1 = pp.omega(end,:);
    end
    
    om.restartMessage();
    
else
    
    time = 0;
    iT = 0;
end

% begin time loop
while time < prams.T
    
    tSingleStep = tic;   

    
    geom = capsules(prams, xc, tau);
    
    [densityF,densityW,Up,wp,stokes,rot,iter,flag,res] = tt.timeStep(geom, tau, ...
                            walls, options, prams);
    
    % update centres and angles
    if (iT == 0 || options.tstep_order == 1) % use forward Euler
        
        xc = xc + tt.dt*Up;
        tau = tau + tt.dt*wp;
        
        if (options.tstep_order == 2)
            Up_m1 = Up;
            wp_m1 = wp;
        end
        
    else if (options.tstep_order == 2) % use Adams Bashforth
            xc = xc + (3/2)*tt.dt*Up - (1/2)*tt.dt*Up_m1;
            tau = tau + (3/2)*tt.dt*wp - (1/2)*tt.dt*wp_m1;
            
            Up_m1 = Up;
            wp_m1 = wp;            
        end
    end
    
    % check for collisions
    [near,~] = geom.getZone(geom,1);
    icollision = geom.collision(near,options.ifmm, options.inear, potF, om);
    
    if (icollision)
        om.writeMessage('WARNING: COLLISION DETECTED');
    end
    
    X = geom.getXY();
    
    time = time + tt.dt;
    iT = iT + 1;
    
    om.writeMessage(....
        ['Finished t=', num2str(time, '%3.3e'), ' in ' num2str(iter) ...
         ' iterations after ', num2str(toc(tSingleStep), '%2i'), ' seconds (residual ', ...
         num2str(res), ')']);
    
    if flag ~= 0
       om.writeMessage(['WARNING: GMRES flag = ', num2str(flag)]); 
    end
        
    om.writeData(time, xc, tau, Up, wp, stokes, rot, densityF, densityW);   
end

Xfinal = X;

om.writeStars();
om.writeMessage(....
        ['Finished entire simulation in ', num2str(toc(ttotal)), ' seconds']);
if (om.profile)
   profsave(profile('info'), om.profileFile);
end


end %rigid2D
