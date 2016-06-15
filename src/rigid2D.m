function [Xfinal] = rigid2D(options, prams, xc, tau)

ttotal = tic;

om = monitor(options, prams);
tt = tstep(options, prams, om);
op = poten(prams.N);

if (om.profile)
    profile on;
end

% read in existing file if necessary
if options.append
    
    pp = post(['../output/data/',options.fileBase, '.dat'], ...
    ['../output/data/',options.fileBase,'_density.dat']);

    time = pp.times(end);
    iT = length(pp.times);
    
    xc = [pp.centres_x(end,:); pp.centres_y(end,:)];
    tau = pp.orientations(end,:);
    
    if (options.tstep_order == 2)
        Up_m1 = [pp.u_x(end,:); pp.u_y(end,:)];
        wp_m1 = pp.omega(end,:);
    end
    
    om.restartMessage();
    
else
    om.writeData(0, xc, tau, zeros(1,prams.nv), zeros(1,prams.nv), zeros(1,prams.nv));
    time = 0;
    iT = 0;
end

% begin time loop
while time < prams.T
    
    tSingleStep = tic;
    
    time = time + tt.dt;
    iT = iT + 1;
    
    geom = capsules(prams, xc, tau);
    
    [density,Up,wp,iter,flag, res] = tt.timeStep(geom);
    
    % update centres and angles
    if (iT == 1 || options.tstep_order == 1) % use forward Euler
        
        xc = xc + tt.dt*Up;
        tau = tau + tt.dt*wp;
        
        if (options.tstep_order == 2)
            Up_m1 =  Up;
            wp_m1 = wp;
        end
        
    else if (options.tstep_order == 2) % use Adams Bashforth
            xc = xc + (3/2)*tt.dt*Up - (1/2)*tt.dt*Up_m1;
            tau = tau + (3/2)*tt.dt*wp - (1/2)*tt.dt*wp_m1;
            
            Up_m1 =  Up;
            wp_m1 = wp;            
        end
    end
    
    [near,~] = geom.getZone(geom,1);
    icollision = geom.collision(near,options.ifmm, options.inear, op, om);
    
    if (icollision)
        om.writeMessage('WARNING: COLLISION DETECTED');
    end
    
    X = geom.getXY();

    om.writeMessage(....
        ['Finished t=', num2str(time, '%3.3e'), ' in ' num2str(iter) ...
         ' iterations after ', num2str(toc(tSingleStep), '%2i'), ' seconds (residual ', ...
         num2str(res), ')']);
    
    if flag ~= 0
       om.writeMessage(['WARNING: GMRES flag = ', num2str(flag)]); 
    end
    
    om.writeData(time, xc, tau, Up(1,:), Up(2,:), wp);
    om.writeDensity(time, density);      
end

Xfinal = X;

om.writeStars();
om.writeMessage(....
        ['Finished entire simulation in ', num2str(toc(ttotal)), ' seconds']);
if (om.profile)
   profsave(profile('info'), om.profileFile);
end

end %rigid2D
