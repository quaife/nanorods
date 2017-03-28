function [Xfinal] = rigid2DCollisions(options, prams, xc, tau, xWalls)

ttotal = tic;
    
geom = capsules(prams, xc, tau);

if (options.confined)
    walls = capsules(prams, xWalls);
else
    walls = [];
    prams.Nw = 0;
    prams.nw = 0;
end

om = monitor(options, prams, xc, tau);
tt = tstep(options, prams, om, geom, walls, tau);
potP = poten(geom.N,om);

if (om.profile)
    profile on;
end

time = 0;
iT = 0;

% begin time loop
while time < prams.T
    
    tSingleStep = tic;   
    geom = capsules(prams, xc, tau);
    X = geom.getXY();
    
    [xc,tau,densityF,densityW,Up,wp,stokes, rot, ...
        ~,~,iter,flag,res] = tt.timeStep(geom, walls, xc, tau);

    time = time + tt.dt;
    iT = iT + 1;
    
    % check for collisions
    [near,~] = geom.getZone(geom,1);
    icollision = geom.collision(near,options.fmm, options.near_singular, potP, om);
    
    if (icollision)
        om.writeMessage('WARNING: COLLISION DETECTED');
    end
    
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
