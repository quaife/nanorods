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

% set minSep
if prams.np > 0
  oc = curve;
  [~,len] = oc.geomProp(geom.X);
  prams.minimum_separation = prams.minimum_separation*max(len)/prams.Np;
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
    
    if options.display_solution
        fill(geom.X(1:end/2,:),geom.X(end/2+1:end,:),'k');

        axis equal
%         xlim([-2.5,2.5])
%         ylim([-2.5,2.5])
        
        if options.confined
            hold on
            for k = 1:prams.nw
                plot(xWalls(1:end/2,k),xWalls(end/2+1:end,k), 'r', 'linewidth', 2);
            end
            
            hold off
        end
        drawnow
    end
    
    if (iT == 0)
        [xc,tau,densityF,densityW,Up,wp,stokes,rot, ...
                ~,~,iter,flag,res] = tt.timeStep(geom, walls, xc, tau, [], [], true);
    else
        [xc,tau,densityF,densityW,Up,wp,stokes,rot, ...
                ~,~,iter,flag,res] = tt.timeStep(geom, walls, xc, tau, Up, wp, false);
    end

    time = time + tt.dt;
    iT = iT + 1;
    
%     % check for collisions
%     [near,~] = geom.getZone(geom,1);
%     icollision = geom.collision(near,options.fmm, options.near_singular, potP, om);
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

Xfinal = X;

om.writeStars();
om.writeMessage(....
        ['Finished entire simulation in ', num2str(toc(ttotal)), ' seconds']);
if (om.profile)
   profsave(profile('info'), om.profileFile);
end

disp(['Total matvecs:', num2str(tt.matvecs)]);

end %rigid2D
