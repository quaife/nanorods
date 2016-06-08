function [Xfinal] = rigid2D(options, prams, xc, tau)

ttotal = tic;

om = monitor(options, prams);
tt = tstep(options, prams, om);

if (om.profile)
    profile on;
end

om.writeData(0, xc, tau, zeros(1,prams.nv), zeros(1,prams.nv), zeros(1,prams.nv));

time = 0;
iT = 0;

while time < prams.T
    
    tic;
    
    time = time + tt.dt;
    iT = iT + 1;
    
    geom = capsules(prams, xc, tau);
    
    [density,Up,wp,iter,flag, res] = tt.timeStep(geom);
    
    % update centres and angles
    if (options.tstep_order == 2)
        xc_m2 =  xc;
        tau_m2 = tau;
    end
    
    if (iT == 1 || options.tstep_order == 1) % use forward Euler
        
        xc = xc + tt.dt*Up;   
        tau = tau + tt.dt*wp;         
        
    else if (options.tstep_order == 2) % use a second oder explicit BDF
            xc = (4/3)*xc - (1/3)*xc_m2 + (2/3)*tt.dt*Up;
            tau = (4/3)*tau - (1/3)*tau_m2 + (2/3)*tt.dt*wp;
        end
    end
    
    X = geom.getXY();

    om.writeMessage(....
        ['Finished t=', num2str(time), ' in ' num2str(iter) ...
         ' iterations after ', num2str(toc), ' seconds (residual ', ...
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

