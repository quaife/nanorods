function [Xfinal, matvecs] = rigid2DCollisions(options, prams, xc, tau, xWalls, torqueP0, Xparticles)

ttotal = tic;
om = monitor(options, prams, xc, tau);

if size(tau,1) > 1
    restart = true;
else
    restart = false;
end

if ~isempty(xc)
    xc = xc(:,:,end);
    tau = tau(end,:);
else
    xc = [];
    tau = [];
end


if nargin == 7
    geom = capsules(prams, Xparticles);
else
    geom = capsules(prams, xc, tau);
end

%prams.minimum_separation = 0.05;

if restart
    min_sep = prams.minimum_separation;
    Np = prams.Np;
    Nw = prams.Nw;
    
    load(['../output/data/', options.file_base, '.mat']);

    end_index = length(t)-1;
    time = t(end_index);
    iT = end_index;
    prams.minimum_separation = min_sep;
    prams.Np = Np;
    prams.Nw = Nw;
%     prams.T = 20;
    %prams.number_steps = 200;
    
    if time > 0
        xc = xc(:,:,end_index);
        tau = tau(end_index,:);
        Up = U(:,:,end_index);
        wp = omega(end_index,:);
        forceP = force_p(:,:,end_index);
        torqueP = torque_p(:,end_index);
        densityF = eta_p(:,:,end_index);
        
        geomOld = capsules(prams, xc,tau);
        
        if options.confined
            densityW = eta_w(:,:,end_index);
        else
            densityW = [];
        end
    end
else
    time = 0;
    iT = 0;
end

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
    
   %sep_w = 0;
   prams.minimum_separation = prams.minimum_separation*max(sep_p, sep_w);
end

%prams.minimum_separation =  0.0320;
%prams.minimum_separation = 0.0638;
%prams.minimum_separation = 3.131636582621994e-03;
options.sedementation = false;
tt = tstep(options, prams, om, geom, walls, tau, torqueP0);

if (om.profile)
    profile on;
end

step_number = 1;
dt0 = tt.dt;
disp(['Initial time step: ', num2str(tt.dt)]);

% begin time loop
while time < prams.T
    
    if mod(step_number,50) == 0 && tt.dt < dt0
        tt.dt = dt0;
        om.writeMessage(['Resetting dt to: ', num2str(tt.dt)]);
    end
    
    tSingleStep = tic;
    geom = capsules(prams, xc, tau);
    X = geom.getXY();
    
    %options.display_solution =true;
    if options.display_solution
        subplot(2,2,1)
        fill(geom.X(1:end/2,:),geom.X(end/2+1:end,:),'k');

        if options.confined
            xmin = min(min(xWalls(1:end/2,:)));
            xmax = max(max(xWalls(1:end/2,:)));
            ymin = min(min(xWalls(end/2+1:end,:)));
            ymax = max(max(xWalls(end/2+1:end,:)));
        else
            xmin = min(min(geom.X(1:end/2,:)));
            xmax = max(max(geom.X(1:end/2,:)));
            ymin = min(min(geom.X(end/2+1:end,:)));
            ymax = max(max(geom.X(end/2+1:end,:)));
        end
        
        axis equal
%         xmin = -5;
%         xmax = 5;
%         ymin = -5;
%         ymax = 5;
        
        xlim([xmin, xmax])
        ylim([ymin, ymax])
        
        title(['t = ',num2str(time)]);
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
                forceP,torqueP,iter,flag,res] = tt.timeStep(geom, [], walls, ...
                xc, tau, [], [], true, forceP, torqueP, densityF, densityW, time);
    else
        [xc,tau,densityF,densityW,Up,wp,stokes,rot, ...
            forceP,torqueP,iter,flag,res] = tt.timeStep(geom, geomOld, walls, ...
            xc, tau, Up, wp, false, forceP, torqueP, densityF, densityW, time);
    end
    
    geomOld = geom;
    
    time = time + tt.dt;
    iT = iT + 1;
    
    % check for collisions
    %if ~options.resolve_collisions
    %    [near,~] = geom.getZone(geom,1);
    %    icollision = geom.collision(near,options.fmm, options.near_singular, potP, om);
        
    %    if (icollision)
    %        om.writeMessage('WARNING: COLLISION DETECTED');
    %    end
    %end
    
    om.writeMessage(....
        ['Finished t=', num2str(time, '%3.3e'), ' in ' num2str(iter) ...
         ' iterations after ', num2str(toc(tSingleStep), '%2i'), ' seconds (residual ', ...
         num2str(res), ')']);
    
    if flag ~= 0
       om.writeMessage(['WARNING: GMRES flag = ', num2str(flag)]); 
    end
        
    om.writeData(time, xc, tau, Up, wp, stokes, rot, densityF, densityW, forceP, torqueP);   
    
    step_number = step_number + 1;
end

Xfinal = X;

om.writeStars();
om.writeMessage(....
        ['Finished entire simulation in ', num2str(toc(ttotal)), ' seconds']);
if (om.profile)
   profsave(profile('info'), om.profileFile);
end

disp(['Total matvecs:', num2str(tt.matvecs)]);
matvecs = tt.matvecs;

end %rigid2D
