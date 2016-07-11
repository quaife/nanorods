function [Xfinal] = rigid2DPeriodic(options, prams, xc, tau, xmin, xmax, ymin, ymax, buffer)

% delete(gcp('nocreate'))
% parpool(options.n_cores_matlab)   

ttotal = tic;
    
geom = capsules(prams, xc, tau);
om = monitor(options, prams);
op = poten(om);

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
    wp = pp.omega(end,:);
    
    if (options.tstep_order == 2)
        Up_m1 = [pp.u_x(end,:); pp.u_y(end,:)];
        wp_m1 = pp.omega(end,:);
    end
    
    om.restartMessage();
    
else
    om.writeData(0, xc, tau, zeros(1,prams.nv), zeros(1,prams.nv), zeros(1,prams.nv));
    time = 0;
    iT = 0;
        
    wp = zeros(size(tau));
end

% begin time loop
while time < prams.T
    
    tSingleStep = tic;

    iT = iT + 1;
    
    tt = tstep(options, prams, om, geom);
    
    time = time + tt.dt;
    [density,Up,wp,iter,flag, res] = tt.timeStep(geom, tt.dt*wp);
    
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
    
    % apply periodic boundary conditions
    X = geom.getXY();
    max_attempts = 50;
    
    for i = 1:prams.nv
        
        % move fibre
        if (xc(1,i) > xmax + buffer && xc(2,i) > 0)
            
            [xc, at] = enforce_periodicity(xmin - max(prams.lengths),...
                xmin, 0, ymax, prams, xc, tau, i, max_attempts);
            
            if at ~= max_attempts
                om.writeMessage(['PERIODIC BOUNDARY CONDITIONS SUCCESSFULLY ENFORCED AFTER ', num2str(at), ' ATTEMPTS']);
                
            else 
                om.writeMessage('PERIODIC BOUNDARY CONDTIONS NOT ENFORCED, TRYING LARGER RANGE');
                [xc, at] = enforce_periodicity(xmin - 2*max(prams.lengths),...
                    xmin, 0, ymax, prams, xc, tau, i, max_attempts);
                
                if at ~= max_attempts
                    om.writeMessage(['PERIODIC BOUNDARY CONDITIONS SUCCESSFULLY ENFORCED AFTER ', num2str(at), ' ATTEMPTS']);
                else
                   om.writeMessage('UNABLE TO ENFORCE PERIODIC BOUNDARY CONDITIONS'); 
                end
            end
            
        else if (xc(1,i) < xmin - buffer && xc(2,i) < 0)
                
                [xc, at] = enforce_periodicity(xmax,...
                    xmax + max(prams.lengths), ymin, 0, prams, xc, tau, i, max_attempts);
                
                if at ~= max_attempts
                    om.writeMessage(['PERIODIC BOUNDARY CONDITIONS SUCCESSFULLY ENFORCED AFTER ', num2str(at), ' ATTEMPTS']);
                    
                else
                    om.writeMessage('PERIODIC BOUNDARY CONDTIONS NOT ENFORCED, TRYING LARGER RANGE');
                    [xc, at] = enforce_periodicity(xmax,...
                        xmax + 2*max(prams.lengths), ymin, 0, prams, xc, tau, i, max_attempts);
                    if at ~= max_attempts
                        om.writeMessage(['PERIODIC BOUNDARY CONDITIONS SUCCESSFULLY ENFORCED AFTER ', num2str(at), ' ATTEMPTS']);
                    else
                        om.writeMessage('UNABLE TO ENFORCE PERIODIC BOUNDARY CONDITIONS');
                    end
                    
                end
            end
        end
        
            
            %xc(1,i) = xmin - buffer;
           
%             for at = 1:max_attempts 
%                 xc(1,i) = xmin - rand(1)*max(prams.lengths);
%                 xc(2,i) = ymax*rand(1);
%                 
%                 % check for collisions
%                 
%                 geom = capsules(prams, xc, tau);
%                 
%                 [near,~] = geom.getZone(geom,1);
%                 
%                 %icollision = geom.collision(near,options.ifmm, options.inear, op, om);
%                 
%                 %if (icollision)
%                 if (length(near.nearFibers{i} > 0))
%                     om.writeMessage(['WARNING: COLLISION DETECTED AFTER ENFORCING PERIODICITY (attempt ', num2str(at), ')']);
%                     
%                 else
%                     om.writeMessage(['PERIODIC BOUNDARY CONDITIONS SUCCESSFULLY ENFORCED AFTER ', num2str(at) ' ATTEMPTS']);
%                     break;
%                 end
%             end
            
%         else if (xc(1,i) < xmin - buffer && xc(2,i) < 0) 
%                 
%                 %xc(1,i) = xmax + buffer;
%                 
%                 for at = 1:max_attempts
%                     xc(1,i) = xmax + rand(1)*max(prams.lengths);
%                     xc(2,i) = ymin*rand(1);
%                     
%                     % check for collisions
%                     
%                     geom = capsules(prams, xc, tau);
%                     
%                     [near,~] = geom.getZone(geom,1);
%                     %icollision = geom.collision(near,options.ifmm, options.inear, op, om);
%                     
%                     if (length(near.nearFibers{i} > 0))%if (icollision)
%                         om.writeMessage(['WARNING: COLLISION DETECTED AFTER ENFORCING PERIODICITY (attempt ', num2str(at), ')']);
%                         
%                     else
%                         om.writeMessage(['PERIODIC BOUNDARY CONDITIONS SUCCESSFULLY ENFORCED AFTER ', num2str(at) ' ATTEMPTS']);
%                         break;
%                     end
%                 end
%             end
%         end
    end    
    
    % check for collisions after updating positions
    geom = capsules(prams, xc, tau);
    
    [near,~] = geom.getZone(geom,1);
    icollision = geom.collision(near,options.ifmm, options.inear, op, om);
    
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

function [xc,at] = enforce_periodicity(xmin, xmax, ymin, ymax, prams, xc, tau, i, max_attempts)

for at = 1:max_attempts
    xc(1,i) = xmin + rand(1)*(xmax - xmin);
    xc(2,i) = ymin + rand(1)*(ymax - ymin);
    
    % check for collisions
    
    geom = capsules(prams, xc, tau);
    
    [near,~] = geom.getZone(geom,1);
    
    %icollision = geom.collision(near,options.ifmm, options.inear, op, om);
    
    %if (icollision)
    if (length(near.nearFibers{i}) == 0)
        break;
    end
end

end %enforce_periodicity
