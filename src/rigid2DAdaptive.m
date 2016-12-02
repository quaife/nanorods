function [Xfinal] = rigid2DAdaptive(options, prams, xc, tau)

global Up wp density

ttotal = tic;

om = monitor(options, prams);
tt = tstep(options, prams, om);

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
    om.writeData(0, xc, tau, zeros(1,prams.nv), zeros(1,prams.nv), zeros(1,prams.nv));
    time = 0;
end



options = odeset('RelTol', 1e-5, 'AbsTol', 1e-8, 'InitialStep', tt.dt, 'Stats', 'on');
Xin = [xc(1,:)';xc(2,:)';tau'];

while (time < prams.T)
    
	tLoop = tic;
    
    % call adaptive ODE solver
    [tSteps, Xout] = ode45(@(t,y) compute_velocities(t, y, prams, tt, om),...
                                        [time, time + tt.dt], Xin, options);
    
    Xin = Xout(end,:);
    xc = [Xin(1:end/3);Xin(end/3+1:2*end/3)];
    tau = Xin(2*end/3+1:end);
    
    time = time + tt.dt;
    
    om.writeMessage(....
        ['Finished full time step for t =', num2str(time, '%3.3e'), ' in ',...
                num2str(toc(tLoop)), ' seconds']);
         
    om.writeMessage(['ODE solver used substeps: ', num2str(tSteps')]);

    om.writeData(time, xc, tau, Up(1,:), Up(2,:), wp);
    om.writeDensity(time, density);   
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
function xdot = compute_velocities(t, y, prams, tt, om)

global Up wp density


tic;

xtmp = y(1:end/3);
ytmp = y(end/3 + 1 : 2*end/3);

xc = [xtmp';ytmp'];
tau = y(2*end/3 + 1: end)';

geom = capsules(prams, xc, tau);
[density,Up,wp,iter,flag, res] = tt.timeStep(geom);

xdot = [Up(1,:)'; Up(2,:)'; wp'];

om.writeMessage(....
    ['Finished substep at t=', num2str(t, '%3.3e'), ' in ' num2str(iter) ...
    ' iterations after ', num2str(toc, '%2i'), ' seconds (residual ', ...
    num2str(res), ')']);

if flag ~= 0
    om.writeMessage(['WARNING: GMRES flag = ', num2str(flag)]);
end

end
