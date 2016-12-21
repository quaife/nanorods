function [Xfinal] = rigid2DAdaptive(options, prams, xc, tau, xWalls)

ttotal = tic;

geom = capsules(prams, xc, tau);

if (options.confined)
    walls = capsules(prams, xWalls);
else
    walls = [];
end

om = monitor(options, prams, xc, tau);
tt = tstep(options, prams, om, geom, walls, tau);
potF = poten(geom.N,om);

if (om.profile)
    profile on;
end

% read in existing file if necessary
if options.append
    
    pp = post(['../output/data/',options.fileBase, '.dat'], ...
        ['../output/data/',options.fileBase,'_density.dat']);
    
    time = pp.times(end);
    xc = [pp.centres_x(end,:); pp.centres_y(end,:)];
    tau = pp.orientations(end,:);
    
    om.restartMessage();
    
else
    
    time = 0;
end

while (time < prams.T)
    
    tSingleStep = tic;  
    
    [xc, tau, dt, num_steps, rel_err, densityF, densityW, ...
                stokes, rot, Up, wp, iter, flag, res] = adaptive_rk(prams.a, prams.b1, prams.b2, ...
                                xc, tau, walls, tt, options, prams);
   
    
    geom = capsules(prams, xc, tau);
    
    tt = tstep(options, prams, om, geom, walls, tau);   
    tt.dt = dt; 
    time = time + tt.dt;
    
    om.writeMessage(....
        ['Finished time step for t =', num2str(time, '%3.3e'), ' in ',...
        num2str(toc(tSingleStep)), ' seconds']);
    
    om.writeMessage(['Adaptive RK solver: res = ',num2str(rel_err), ...
        ' , number of step halves = ', num2str(num_steps - 1)]);
    
%      % check for collisions
%     [near,~] = geom.getZone(geom,1);
%     icollision = geom.collision(near,options.ifmm, options.inear, potF, om);
%     
%     if (icollision)
%         om.writeMessage('WARNING: COLLISION DETECTED');
%     end
    
    om.writeMessage(....
        ['Finished t=', num2str(time, '%3.3e'), ' in ' num2str(iter) ...
         ' iterations after ', num2str(toc(tSingleStep), '%2i'), ' seconds (residual ', ...
         num2str(res), ')']);
    
    if flag ~= 0
       om.writeMessage(['WARNING: GMRES flag = ', num2str(flag)]); 
    end
        
    om.writeData(time, xc, tau, Up, wp, stokes, rot, densityF, densityW); 
end

geom = capsules(prams, xc, tau);
Xfinal = geom.getXY();

om.writeStars();
om.writeMessage(....
    ['Finished entire simulation in ', num2str(toc(ttotal)), ' seconds']);
if (om.profile)
    profsave(profile('info'), om.profileFile);
end

end %rigid2D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xc_new, tau_new, dt, num_steps, rel_err, densityF, densityW, ...
                stokes, rot, Up, wp, iter, flag, res] = adaptive_rk(a, b1, b2, ...
                                xc, tau, walls, tt, options, prams)

dt = tt.dt;
rel_err = 2*prams.tol;
num_steps = 0;

while (rel_err > prams.tol && dt > prams.dt_min)
    % compute lower order solution using b1
    
    k = zeros(length(xc(:))+length(tau),length(b1));
    
    for s = 1:length(b1)
        
        dTmp = k(:,1:s-1)*a(s,1:s-1)';
        
        dx = reshape(dTmp(1:2*prams.nv), prams.nv, 2);
        dtau = dTmp(2*prams.nv+1:end);
        
        xc_k = xc + dt*dx';
        tau_k = tau + dt*dtau';
        k(:,s) = compute_velocities(xc_k, tau_k, walls, ...
                                    tt, options, prams);
                                
        if s == 1
            Up = [k(1:prams.nv,s)';k(prams.nv+1:2*prams.nv,s)'];
            wp = k(2*prams.nv+1:end,s)';
        end
    end
    
    sol_low = [xc(1,:)';xc(2,:)';tau'] + dt*k*b1;
    
    % add in last stage to compute higher order solution
    dTmp = k*a(end,:)';
        
    dx = reshape(dTmp(1:2*prams.nv), prams.nv, 2);
    dtau = dTmp(2*prams.nv+1:end);

    xc_k = xc + dt*dx';
    tau_k = tau + dt*dtau';
        
    [k(:,end+1),densityF, densityW, stokes, rot, iter, flag, res] = ...
                compute_velocities(xc_k, tau_k, walls, tt, options, prams);
    
    sol_high = [xc(1,:)';xc(2,:)';tau'] + dt*k*b2;
    
    % compute difference, adjust step size if necessary
    
    rel_err = max(abs(sol_high - sol_low)./sol_high);
    num_steps = num_steps + 1;
    
    if options.verbose
        disp(['Relative error = ', num2str(rel_err)]);
    end
    
    if (rel_err > prams.tol)
        dt = dt/2;
        if options.verbose
        	disp(['Halving time step size, new dt = ', num2str(dt)]);
        end
    end
end

xc_new = [sol_high(1:prams.nv)';sol_high(prams.nv+1:2*prams.nv)'];
tau_new = sol_high(2*prams.nv+1:end)';


end %adaptive_rk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xdot, densityF, densityW, stokes, rot, iter, flag, res] = ...
                    compute_velocities(xc, tau, walls, tt, options, prams)
    
    geom = capsules(prams, xc, tau);
    [densityF,densityW,Up,wp,stokes,rot,iter,flag,res] = tt.timeStep(geom, tau, ...
                            walls, options, prams);
     
     xdot = [Up(1,:)';Up(2,:)';wp'];
end
