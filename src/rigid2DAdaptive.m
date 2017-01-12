function [Xfinal] = rigid2DAdaptive(options, prams, xc, tau, xWalls)

ttotal = tic;
geom = capsules(prams, xc, tau);


if (options.confined)
    walls = capsules(prams, xWalls);
else
    walls = [];
    prams.Nbd = 0;
    prams.nbd = 0;
end

om = monitor(options, prams, xc, tau);
tt = tstep(options, prams, om, geom, walls, tau);
potF = poten(geom.N,om);

if (om.profile)
    profile on;
end

% set Butcher tableau coefficients
switch options.tstep_order
    case 1 % Heun-Euler
        a = [0, 0; 1, 0];
        b1 = 1;
        b2 = [1/2; 1/2];
        
    case 2 % Bogacki-Shampine
        a = [0,0,0,0;
            1/2,0,0,0;
            0,3/4,0,0;
            2/9,1/3,4/9,0];
        b1 = [2/9;1/3;4/9];
        b2 = [7/24;1/4;1/3;1/8];
        
    case 4 % Runge-Kutta-Fehlberg
        a = [0, 0, 0, 0, 0, 0;
             1/4 ,0, 0, 0, 0, 0;
             3/32, 9/32, 0, 0, 0, 0;
             1932/21971, -7200/2197, 7296/2197, 0, 0, 0;
             439/216, -8, 3680/513, -845/4104, 0, 0;
             -8/27, 2, -3544/2565, 1859/4104, -11/40, 0];
         b1 = [25/216; 0; 1408/2565; 2197/4104; -1/5];
         b2 = [16/135; 0; 6656/12825; 28561/56430; -9/50; 2/55];
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
    
    [xc_new, tau_new, accept, rel_err, densityF, densityW, ...
                stokes, rot, Up, wp, iter, flag, res] = adaptive_rk(a, ...
                                b1, b2, xc, tau, walls, tt, potF, om,...
                                options, prams);
    if ~accept
        % halve time step size, try again until accept
        while (~accept && tt.dt > options.dt_min)
            tt.dt = tt.dt/2;
            
            if options.verbose
                om.writeMessage(['Halving dt, new dt = ', num2str(tt.dt)]);
            end
            
            [xc_new, tau_new, accept, rel_err, densityF, densityW, ...
                stokes, rot, Up, wp, iter, flag, res] = adaptive_rk(a, ...
                                b1, b2, xc, tau, walls, tt, potF, om,...
                                options, prams);
        end
    end
    
    time = time + tt.dt;
    xc = xc_new;
    tau = tau_new;
    
    % adjust time step size    
    tt.dt =  max(min(tt.dt*(options.rk_tol/rel_err)^(1/(options.tstep_order+1)), ...
                        2*tt.dt), 0.25*tt.dt);                    
    if options.verbose
         om.writeMessage(['Step completed, new dt = ', num2str(tt.dt)]);
    end
    
    % print time step information
    om.writeMessage(....
        ['Finished time step for t =', num2str(time, '%3.3e'), ' in ',...
        num2str(toc(tSingleStep)), ' seconds']);
    
    om.writeMessage(....
        ['Finished t=', num2str(time, '%3.3e'), ' in ' num2str(iter) ...
         ' iterations after ', num2str(toc(tSingleStep), '%2i'), ...
         ' seconds (residual ', num2str(res), ')']);
    
    if flag ~= 0
       om.writeMessage(['WARNING: GMRES flag = ', num2str(flag)]); 
    end
        
    % save data
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
function [xc_new, tau_new, accept, rel_err, densityF, densityW, ...
    stokes, rot, Up, wp, iter, flag, res] = adaptive_rk(a, b1, ...
    b2, xc, tau, walls, tt, potF, om, options, prams)

dt = tt.dt;
accept = true;
rel_err = 1;

% compute lower order solution using b1
k = zeros(length(xc(:))+length(tau),length(b1));

for s = 1:length(b1)
    
    dTmp = k(:,1:s-1)*a(s,1:s-1)';
    
    dx = reshape(dTmp(1:2*prams.nv), prams.nv, 2);
    dtau = dTmp(2*prams.nv+1:end);
    
    xc_k = xc + dt*dx';
    tau_k = tau + dt*dtau';
    
    [k(:,s),densityF, densityW, stokes, rot, iter, flag, res] = ...
                compute_velocities(xc_k, tau_k, ...
                                    walls, tt, options, prams);
                                                                                      
end

Up = [k(1:prams.nv,1)';k(prams.nv+1:2*prams.nv,1)'];
wp = k(2*prams.nv+1:end,1)';    
sol_low = [xc(1,:)';xc(2,:)';tau'] + dt*k*b1;

% check for collisions
geom = capsules(prams, [sol_low(1:prams.nv),sol_low(prams.nv+1:2*prams.nv)]',...
                            sol_low(2*prams.nv+1:end));                   
[near,~] = geom.getZone(geom,1);                       
icollision = geom.collision(near,options.ifmm, options.inear, potF, om);

if (icollision)
    om.writeMessage('WARNING: COLLISION DETECTED USING LOWER ORDER RK');
    accept = false;
end

if accept
    % add in last stage to compute higher order solution
    dTmp = k*a(end,1:end-1)';
    
    dx = reshape(dTmp(1:2*prams.nv), prams.nv, 2);
    dtau = dTmp(2*prams.nv+1:end);
    
    xc_k = xc + dt*dx';
    tau_k = tau + dt*dtau';
    
    [k(:,end+1),densityF, densityW, stokes, rot, iter, flag, res] = ...
        compute_velocities(xc_k, tau_k, ...
                        walls, tt, options, prams);
    
    sol_high = [xc(1,:)';xc(2,:)';tau'] + dt*k*b2;
    
    % check for collisions
    geom = capsules(prams, [sol_high(1:prams.nv),sol_high(prams.nv+1:2*prams.nv)]',...
                            sol_high(2*prams.nv+1:end));
    [near,~] = geom.getZone(geom,1);  
    icollision = geom.collision(near,options.ifmm, options.inear, potF, om);
    
    if (icollision)
        om.writeMessage('WARNING: COLLISION DETECTED USING HIGHER ORDER RK');
        accept = false;
    else

        % compute difference, accept/reject
        err = abs(sol_high - sol_low);
        err(err < 1e-14) = 0;
        
        rel_err = max(err./abs(sol_high));
        
        if options.verbose
            om.writeMessage(['Relative error = ', num2str(rel_err)]);
        end
        
        if rel_err > options.rk_tol
            accept = false;
        end
    end

else
    sol_high = zeros(3*prams.nv,1);
end

xc_new = [sol_high(1:prams.nv)';sol_high(prams.nv+1:2*prams.nv)'];
tau_new = sol_high(2*prams.nv+1:end)';


end %adaptive_rk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xdot, densityF, densityW, stokes, rot, iter, flag, res] = ...
                    compute_velocities(xc, tau, walls, tt, options, prams)
    
    geom = capsules(prams, xc, tau);
    [densityF,densityW,Up,wp,stokes,rot,iter,flag,res] = tt.timeStep(geom, ...
                            tau, walls, options, prams);
     
     xdot = [Up(1,:)';Up(2,:)';wp'];
end
