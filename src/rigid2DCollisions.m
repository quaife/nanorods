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
    sep_p = max(len)/prams.Np;
    
    if ~isempty(xWalls)
        [~,len] = oc.geomProp(walls.X);
        sep_w = max(len)/prams.Nw;
    else
        sep_w = sep_p;
    end
    
    prams.minimum_separation = prams.minimum_separation*max(sep_p, sep_w);
end

%prams.minimum_separation = 0.019991;

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
        subplot(2,2,1)
        fill(geom.X(1:end/2,:),geom.X(end/2+1:end,:),'k');

        axis equal
        xlim([-10,10])
        ylim([-10,10])
        
        if options.confined
            hold on
            for k = 1:prams.nw
                plot(xWalls(1:end/2,k),xWalls(end/2+1:end,k), 'r', 'linewidth', 2);
            end
            
            hold off
        end
        
        if iT > 1
            subplot(2,2,2);
            plot(densityF);
            
            subplot(2,2,3);
            plot(forceP)
            
            subplot(2,2,4);
            plot(torqueP);
        else
            
            subplot(2,2,3);
            plot(zeros(2,prams.np))
            
            subplot(2,2,4);
            plot(zeros(1,prams.np));
        end
        drawnow
    end
    
    if (iT == 0)
        densityF = zeros(2*prams.Np,prams.np);
        densityW = zeros(2*prams.Nw,prams.nw);
        forceP = zeros(2,prams.np);
        torqueP = zeros(1,prams.np);
        
        [xc,tau,densityF,densityW,Up,wp,stokes,rot, ...
                forceP,torqueP,iter,flag,res] = tt.timeStep(geom, walls, xc, tau, [], [], true, forceP, torqueP, densityF, densityW);
    else
        [xc,tau,densityF,densityW,Up,wp,stokes,rot, ...
                forceP,torqueP,iter,flag,res] = tt.timeStep(geom, walls, xc, tau, Up, wp, false, forceP, torqueP, densityF, densityW);
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
        
    om.writeData(time, xc, tau, Up, wp, stokes, rot, densityF, densityW, forceP, torqueP);   
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
